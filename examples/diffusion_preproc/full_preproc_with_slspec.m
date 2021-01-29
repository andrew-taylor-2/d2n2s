%example preproc script for the following case: one folder with topup b0s, one folder
%   with dwis. 
% 
%unlike the other preproc example, this one uses mrtrix for degibbs and
%   rician correction. I'm not sure which is faster.
%
%slspec correction within eddy can be performed using this script
%
%prior to this script, dicoms should be sorted and converted with dcm2niix
%
% dicom_sort('C:\re\test\AD_study','')
% eval(sprintf('!dcm2niix %s',folder))
clear

topup_dir='/Users/andrew/re/test/test_example_script/20180323_084032_Test_ADStudy/dicom/DKI_30_Dir_TOP_UP_PA_23';%in some cases, including '/' at the end can screw this up %spm_select(Inf,'dir','choose fbi dir');
dki_dir='/Users/andrew/re/test/test_example_script/20180323_084032_Test_ADStudy/dicom/DKI_30_Dir_17';%spm_select(Inf,'dir','choose dki dir');

%% make niftis -- "nii-ify"
dcm2niix_path='/Applications/MRIcroGL/dcm2niix';
niiify=@(x) system([dcm2niix_path ' ' x]); 

niiify(topup_dir);
niiify(dki_dir);

%% grab nii file
%hope there's only one in each folder
fnify=@(x) [x.folder filesep x.name];

topnii_dir=dir([topup_dir filesep '*.nii']);
topnii=fnify(topnii_dir);

dkinii_dir=dir([dki_dir filesep '*.nii']);
dkinii=fnify(dkinii_dir);

%% perform a couple of the basic preproc steps on both datasets (separately)

%*see bottom for a potentially more intuitive way to do this (if you don't
%   like anonymous functions)

%I like this way better :)   (it actually facilitates easy swapping of order and removing
%   certain commands)
%anonymous function that's going to help with filenames
appendd=@(x,str) [choose_output(@() fileparts(x),1) filesep choose_output(@() fileparts(x),2) str choose_output(@() fileparts(x),3)]; %append str to end of file (but before extension)
newfile=@(x,name) [choose_output(@() fileparts(x),1) filesep name choose_output(@() fileparts(x),3)]; %make new file in same folder as first input
gunziped=@(x) x(1:end-3); %for zipped filename strings; removes '.gz'


field_strs={'top','dki'}; %still doesn't seem like the best way bc these strings aren't necessarily related to the variable names... but whatever...
fn={topnii,dkinii};
for i=1:2 
    [path,name,ext]=fileparts(fn{i});
    
    in=fn{i};
    out=in;
    
    conditionn= i==2 || numel(spm_vol(fn{i}))>1; %don't do this if there's only one volume. also, don't even check if it's dwis, that wastes time
    %but also it's super worth considering whether or not we want to throw
    %   undenoised revPE volumes into topup with corrected dwi data
    
    if conditionn
    %denoise
    out=appendd(in,'_dn');
    noise_out=newfile(in,'noise');
    command=['dwidenoise ' in ' ' out ' -noise ' noise_out ' -force'];
    [a,b]=system(command)
    
    %get res
    noise_residuals_out=newfile(in,'res');
    command=['mrcalc ' in ' ' out ' -subtract ' noise_residuals_out];
    [a,b]=system(command)
    end
    
    %gibbs unring
    in=out;
    out=appendd(in,'_ur');
    command=['mrdegibbs ' in ' ' out ' -force'];
    [a,b]=system(command)
    
    if conditionn
    %rician
    %yes, I realize this is a disgusting way to implement this
    %   but it goes with the form of the above commands for simplicity
    in=out;
    out=appendd(in,'_rc');
    signal_sqr=newfile(in,'_sqr');
    noise_sqr=newfile(in,'noise_sqr');
    diff=[path filesep 'diff' ext];
    
    c1=['fslmaths ' in ' -sqr ' signal_sqr ';'];
    c2=['fslmaths ' noise_out ' -sqr ' noise_sqr ';'];
    c3=['fslmaths ' signal_sqr ' -sub ' noise_sqr ' ' diff ';'];
    c4=['fslmaths ' diff ' -sqrt ' out];
    [a,b]=system([c1 c2 c3 c4]);
    end
    
    %handle gzip
    if exist([out '.gz'],'file')
        gunzip([out '.gz']); 
        delete([out '.gz']);
    end
    
    
    %assign the outname to a descriptively named variable
    out_fn.(field_strs{i})=out;
    
end


%% prepare to stick things into topup/eddy
dki_final=out_fn.dki; %the less programmatic way: appendd(dkinii,'_dn_ur_rc'); 
topup_final=out_fn.top; %appendd(topnii,'_dn_ur_rc');

picks.dwis=dki_final;
picks.topups=topup_final; %the way this is put into topup_n_eddy is sort of the opposite of how d2n2s obj.picks goes 

eddy_intermediate_dir=fullfile(dki_dir,'..','eddy_intermediate');
arbitrary_output_name='eddyd'; %will be assigned such that the final image from eddy is [eddy_intermediate_dir filesep 'eddy' filesep arbitrary_output_name '.nii.gz']
%i made this able to be assignable bc I don't want to force people to use
%	the name 'eddyd' lawl
%also it's more transparent in the main body of this script exactly what
%   the "final output image" is 
topup_n_eddy2(dki_dir,topup_dir,eddy_intermediate_dir,picks,arbitrary_output_name); %I'm using a modified version of topup_n_eddy here which includes slspec 
%I could also make a topup_n_eddy that works from objects... enticing...

%coregistration used to go here, but it has been replaced by eddy
%% read data for final organization
%NOTE: d2n2s doesn't really shine here... what I need is a topup_n_eddy which
%   will work with objects.

%change the name of outputs so that d2n2s will work
eddy_final_out_gz=[eddy_intermediate_dir filesep 'eddy' filesep arbitrary_output_name '.nii.gz'];
gunzip(eddy_final_out_gz); %for spm
eddy_final_out=   [eddy_intermediate_dir filesep 'eddy' filesep arbitrary_output_name '.nii'];

eddy_bvecs_out=[eddy_intermediate_dir filesep 'eddy' filesep arbitrary_output_name '.rotated_bvecs'];
copyfile(eddy_bvecs_out,[eddy_bvecs_out '.bvec']) %this is only going to work in general if there are no other .bvecs in the folder

dkibvec_dir=dir([dki_dir filesep '*.bval']); %grab the original dwi bvals, which are unchanged so far
copyfile(fnify(dkibvec_dir),[eddy_intermediate_dir filesep 'eddy' filesep arbitrary_output_name '.bval'])

%use d2n2s
flagg.pick=eddy_final_out;
rdkis=d2n2s(eddy_intermediate_dir,flagg);

%% average b0s:
%find indices of b0s and "read" images
%in all previous tests, using at_read_vols() is the same as just
%   using the .img data, but we're doing this just to be safe
b0i=find([rdkis.bval]==0);
allb0s=at_read_vols(rdkis(b0i)); %given the comment above, why don't I just use cat(4,rdkis(...).img)?
b0sm=mean(allb0s,4);
%we're throwing out corrected topup b0s -- you could also use these, if you
%   want

rdkis(b0i(1)).img=b0sm; %first b0 replaced with mean
rdkis(b0i(2:end))=[]; %other b0s deleted


%% change dke_parameters and write stuff
dke_dir=[dki_dir filesep '..' filesep 'dke'];
mkdir(dke_dir)
flags=[];
d2n2s_write(rdkis,dke_dir,'dkis',flags)

%% analyze
% command = ['C:\Programs\DKE\win64_internal\Executables\dke.exe ' dki_dir '\dke\dke_params.txt'];
% [status,cmdout] = system(command);
% 
% optimize_FBI([fbi_dir '\final'])


%sick function I got from mathworks
function varargout = choose_output(f, out_indices)
% varargout = output(f, out_indices)
% 
% Select outputs from a function with multiple outputs. This is most useful
% when writing anonymous functions, where it is otherwise difficult to use
% anything but the first output.
% 
% Example:
%
% >> index = output(@() min([3 1 4]), 2)
% index =
%      2
%
% Multiple indices work as well.
%
% >> [index, minimum] = output(@() min([3 1 4]), [2 1])
% index =
%      2
% minimum =
%      1
%
% Tucker McClure
% Copyright 2013 The MathWorks, Inc.
    [outs{1:max(out_indices)}] = f();
    varargout = outs(out_indices);
end

function topup_n_eddy2(dwi_dir,topup_dir,outdir,pick,arbitrary_output_name)
% to keep "pick" but remove slspec functionality, simply delete the
% '--slspec=...' in the eddy command below. The earlier code that allows for this
% is not likely to cause errors.
%
% use d2n2s to run eddy from 2 dwi objects
% this assumes you have dwis in a folder and topup b0s acquired differently
% in another folder

%notes: this function will not work if you have spaces in your pathname

%% handle variable input conditions
% THIS SECTION UNDER CONSTRUCTION
%make output dir if it doesn't exist
if ~exist(outdir,'dir')
    mkdir(outdir)
end
% 
% %if something's a nii.gz, convert (to a unique name)
% %check and see if there's .niis in the folders
% originaldir{1}=dwi_dir; %this is used if conversion is needed -- this way of doing it will work for symbolic linked paths
% originaldir{2}=topup_dir;
% checkniis{1}=dir([dwi_dir filesep '*.nii']);
% checkniis{2}=dir([topup_dir filesep '*.nii']);
% checkniigzs{1}=dir([dwi_dir filesep '*.nii.gz']);
% checkniigzs{2}=dir([topup_dir filesep '*.nii.gz']);
% 
% for i=1:2
%     if isempty(checkniis{i}) && isempty(checkniigzs{i})
%         error('images missing -- no .nii or .nii.gz files found for input dir(s)')
%     elseif isempty(checkniis{i}) && ~isempty(checkniigzs{i})
%         if length(checkniigzs{i})>1
%             error('not sure which gz to convert')
%         end
%         gunzip([originaldir{i} filesep checkniigzs{i}.name],);
%         

% if a pick flag exists, use this to choose niftis inputs in a folder with
% multiple niftis


%assign to flags for reading
dflags.pick=[];
tflags.pick=[];

if exist('pick','var') && isfield(pick,'dwis') % could get a list of fieldnames and see if any of the strings contain 'dwi' lol
    dflags.pick=pick.dwis;
end

if exist('pick','var') && isfield(pick,'topups') 
    tflags.pick=pick.topups;
end

%% read data
% dwi_dir='';
% topup_dir='';
dwis=d2n2s(dwi_dir,dflags); %in the future, i might run flags.no='img'. This would necessitate taking dwis.bval==... and putting it into fslroi. It would also allow one to run this script from .nii.gzs
topups=d2n2s(topup_dir,tflags); %note that topup acquisitions might not have accompanying bvals and bvecs like other image series

% read json data
slspec_d_fn=[outdir filesep 'slspec_dwis.txt'];
% Reading topup image slspec data is not necessary bc these images don't become the outputs directly
getslspec(dwis(1).json,slspec_d_fn);

%% make a mask
%find the first b0
idx_b0s=find([dwis.bval]==0);
idx_b0s_1st=idx_b0s(1);

%also do this for topups for later
if isfield(topups,'bval')
    idx_b0s_topups=find([topups.bval]==0);
else %assume that all images are b0s
    idx_b0s_topups=1:numel(topups); %there's probably a better way than to use this; let's see how it shakes out
end

%write to a different folder
b01=dwis(idx_b0s_1st);
intermediate_dir=[outdir filesep 'pre_topup'];
mkdir(intermediate_dir)
d2n2s_write(b01,intermediate_dir,'b0',[])

%use bet
in=[intermediate_dir filesep 'b0.nii'];
mask_fn=[intermediate_dir filesep 'b0_mask.nii'];
command=['bet ' in ' ' in ' -m -n']; % these options will create a brain mask, but no skull stripped b0
[a,b]=system(command); 

%% prepare things for topup proper
%make a topup directory
topup_interm_dir=[outdir filesep 'topup'];
mkdir(topup_interm_dir)

%% use the simple_acqp function to make index and acqp text files
index_fn=[topup_interm_dir filesep 'index.txt'];
acqp_fn=[topup_interm_dir filesep 'acqp.txt'];
simple_acqp(dwis,topups,topup_interm_dir); %this command makes a two line acqp and a corresponding index file

%% grab ALL the b0s for topup
b0s_dwis=dwis(idx_b0s); %use idx from earlier
b0s_topups=topups(idx_b0s_topups);
all_b0s=join_obj(b0s_dwis,b0s_topups);

%write to a topup dir
d2n2s_write(all_b0s,topup_interm_dir,'all_b0s',[]);
topup_b0s_fn=[topup_interm_dir filesep 'all_b0s.nii']; %d2n2s_write makes this file

%since we're not doing any corrections, we don't actually need to write out
%our main 4d image again -- we were just grabbing b0s and metadata. If, for
%some reason, dcm2niix wrote your bvals as tall instead of long, you could
%use d2n2s_write with the .gr = 1 or 0 flag, i forgot which

%but anyway yeah you can grab bvals and bvecs for eddy straight from the
%dcm2niix conversion

%% get bvec, bval, nii files for main 4d image
%pretty sure no other bvec,bval files have been written to this main folder
bvecs=dir([dwi_dir filesep '*.bvec']);
bvecs_fn=[dwi_dir filesep bvecs.name];%hope dir only grabbed one file

bvals=dir([dwi_dir filesep '*.bval']);
bvals_fn=[dwi_dir filesep bvals.name]; %only one plz

if exist('pick','var') && isfield(pick,'dwis') % could get a list of fieldnames and see if any of the strings contain 'dwi' lol
    dwi_image_fn=pick.dwis;
else
    dwi_image_dir=dir([dwi_dir filesep '*.nii']);
    dwi_image_fn=[dwi_dir filesep dwi_image_dir.name]; %only one plz
end


%% run the FSL commands!
eddy_dir=[outdir filesep 'eddy'];
out_fn=[eddy_dir filesep arbitrary_output_name];
mkdir(eddy_dir)


command=['topup --imain=' topup_b0s_fn ' --datain=' acqp_fn ' ' ...
    '--config=b02b0.cnf --out=' topup_interm_dir filesep 'my_topup_results ' ...
    '--fout=' topup_interm_dir filesep 'my_field --iout=' topup_interm_dir filesep ...
    'my_unwarped_images'];
system(command)

command=['eddy --imain=' dwi_image_fn '  --mask=' mask_fn ' --index=' index_fn ...
    ' --acqp=' acqp_fn ' --bvecs=' bvecs_fn ' --bvals=' bvals_fn ' --out=' ...
    out_fn ' --topup=' topup_interm_dir filesep 'my_topup_results --slspec=' slspec_d_fn ' --repol --data_is_shelled'];
system(command)



% function simple_acqp(json_fn,out_dir,num_imgs)
function simple_acqp(dwi_obj,topup_obj,out_dir)
%purpose: build acqp.txt and index.txt

%% find out how many b0s you had from each scan
%this sort of duplicates efforts from the main body of the script...
idx_b0_dwi=find([dwi_obj.bval]==0);
num_dwi_b0s=numel(idx_b0_dwi);

if isfield(topup_obj,'bval')
    idx_b0_topups=find([topup_obj.bval]==0);
else %assume they're all b0s
    idx_b0_topups=1:numel(topup_obj);
end
num_top_b0s=numel(idx_b0_topups);

%% other stuff
% j_struct{1}=readJson(json_fn); %json1);
j_struct{1}=dwi_obj(1).json; %json1);
j_struct{2}=topup_obj(1).json;

%see if you have the fields you need
for i=1:2
    if ~isfield(j_struct{i},'TotalReadoutTime') || ~isfield(j_struct{i},'PhaseEncodingDirection')
        error('can''t find the required fields')
    end
end

%build two lines for acqp.txt
stringss={'j','j-'}; %known direction strings
vectorss={[0 1 0],[0 -1 0]};% corresponding vectors; 

for i=1:2
    time{i}=j_struct{i}.TotalReadoutTime;
    dir_string{i}=j_struct{i}.PhaseEncodingDirection;
    
    if all(dir_string{i}==stringss{1})
        vec{i}=vectorss{1};
    elseif all(dir_string{i}==stringss{2})
        vec{i}=vectorss{2};
    else
        error('Unknown PEdirection string')
    end
    
    line{i}=[vec{i} time{i}]; 
end


% there's one unique line for dwis and one for topups; this is going to
% write them as many times as each set has b0s (bc that's what topup wants)
fidd = fopen([out_dir filesep 'acqp.txt'],'w');
for iii=1:num_dwi_b0s
    fprintf(fidd,'%i %i %i %f\n',line{1}');
end
for jjj=1:num_top_b0s
    fprintf(fidd,'%i %i %i %f\n',line{2}');
end
fclose(fidd);
    


%also make index file
fidd = fopen([out_dir filesep 'index.txt'],'w');

index_array1=repmat(1,[numel(dwi_obj),1]); %this just writes as many 1s as there are dwis
fprintf(fidd,'%i\n',index_array1');

%you might use the next 3 lines if you were to append the topup b0s to the
%end of the dwi set
% index_of_second_unique_acqp_line=num_top_b0s+1;
% index_array2=repmat(index_of_second_unique_acqp_line,[numel(topup_obj),1]); %this just writes as many of the index of the second unique acqp line as there are topups
% fprintf(fidd,'%i\n',index_array2');

fclose(fidd);

end
end

%*different way
% %perform a couple of the basic steps on both datasets (separately)
% for fn={topup_dir,dki_dir} 
%     [path,name,ext]=fileparts(fn{1});
%     
%     %denoise
%     in=fn{1};
%     out=[path filesep name '_dn' ext];
%     noise_out=[path filesep 'noise' ext];
%     command=['dwidenoise ' in ' ' out ' -noise ' noise_out ' -force'];
%     [a,b]=system(command);
%     
%     noise_residuals_out=[path filesep 'res' ext];
%     command=['mrcalc ' in ' ' out ' -subtract ' noise_residuals_out];
%     [a,b]=system(command);
%     
%     %gibbs unring
%     in=[path filesep name '_dn' ext]; % or in=out;
%     out=[path filesep name '_dn_ur' ext];
%     command=['mrdegibbs ' in ' ' out ' -force'];
%     [a,b]=system(command);
%     
%     
%     %rician
%     %yes, I realize this is a disgusting way to implement this
%     %but it goes with the form of the above commands for simplicity
%     in=out;
%     out=[path filesep name '_dn_ur_rc' ext];
%     signal_sqr=[path filesep name '_sqr' ext];
%     noise_sqr=[path filesep 'noise_sqr' ext];
%     diff=[path filesep 'diff' ext];
%     
%     c1=['fslmaths ' in ' -sqr ' signal_sqr ';'];
%     c2=['fslmaths ' path filesep 'noise -sqr ' noise_sqr ';'];
%     c3=['fslmaths '  signal_sqr ' -sub ' noise_sqr ' ' diff ';'];
%     c4=['fslmaths ' diff ' -sqrt ' out];
%     [a,b]=system([c1 c2 c3 c4]);
%     
% end