clear

%% Define Variables 
appendd=@(x,str) [choose_output(@() fileparts(x),1) filesep choose_output(@() fileparts(x),2) str choose_output(@() fileparts(x),3)]; %append str to end of file (but before extension)
newfile=@(x,name) [choose_output(@() fileparts(x),1) filesep name choose_output(@() fileparts(x),3)]; %make new file in same folder as first input
gunziped=@(x) x(1:end-3);
fnify=@(x) [x.folder filesep x.name];

%% grab relevant images and filter 
grab_4ds=dir('/Users/andrew/re/coop/final_run/PD/VaP*/original_dki.nii');

excluded_for_QC=[1,6,8,11,12,14,15,20,22,23,24,28,32];
already_finished=[1,2,3];
dont_do=union(excluded_for_QC,already_finished);
grab_4ds(dont_do)=[];

for i1=1:length(grab_4ds) %already run: 1
tic
try diary([grab_4ds(i1).folder filesep 'diary.txt']);catch;end
    
%% run basic prerproc
in=fnify(grab_4ds(i1));

out=appendd(in,'_dn');
noise_out=newfile(in,'noise');
command=['dwidenoise ' in ' ' out ' -noise ' noise_out ' -force'];
[a,b]=system(command)

%get res
noise_residuals_out=newfile(in,'res');
command=['mrcalc ' in ' ' out ' -subtract ' noise_residuals_out];
[a,b]=system(command)


% no gibbs
% %gibbs unring
% in=out;
% out=appendd(in,'_ur');
% command=['mrdegibbs ' in ' ' out ' -force'];
% [a,b]=system(command)

%rician
%yes, I realize this is a disgusting way to implement this
%   but it goes with the form of the above commands for simplicity
in=out;
out=appendd(in,'_rc');
signal_sqr=newfile(in,'_sqr'); %signal_sqr=appendd(in,'_sqr');
noise_sqr=newfile(in,'noise_sqr');
diff=newfile(in,'diff');

c1=['fslmaths ' in ' -sqr ' signal_sqr ';'];
c2=['fslmaths ' noise_out ' -sqr ' noise_sqr ';'];
c3=['fslmaths ' signal_sqr ' -sub ' noise_sqr ' ' diff ';'];
c4=['fslmaths ' diff ' -sqrt ' out];
[a,b]=system([c1 c2 c3 c4]);

%handle gzip
if exist([out '.gz'],'file')
    gunzip([out '.gz']); 
    delete([out '.gz']);
end

   

flags1.pick=out;
basefolder=choose_output(@() fileparts(out),1);
dkis=d2n2s(basefolder,flags1);

%% eddy bit
%d2n2s to put things into eddy
dwi_dir=[basefolder filesep 'pre_eddy'];
outdir=basefolder; %maybe this will turn out to be a confusing choice, we'll see

d2n2s_write(dkis,dwi_dir,'dkis',[])
just_eddy(dwi_dir,outdir,[])

%grab .nii,.bvec,.bval
eddy_output_folder=[outdir filesep 'eddy']; %these names are based on how "just_eddy" writes things
eddy_final_image=[eddy_output_folder filesep 'eddyd.nii.gz']; %could remove the hard-coding on this sometime
gunzip(eddy_final_image)
eddy_final_image=gunziped(eddy_final_image);

%original bvals
bvals_fn=dir([basefolder filesep '*.bval']); 
copyfile(fnify(bvals_fn),[eddy_output_folder filesep 'eddyd.bval'])

%rename bvec for pickup
movefile([eddy_output_folder filesep 'eddyd.eddy_rotated_bvecs'],[eddy_output_folder filesep 'eddyd.eddy_rotated_bvecs.bvec'])

%d2n2s
ff.pick=eddy_final_image;
corrected_dkis=d2n2s(eddy_output_folder,ff);

%% choosing certain images by their bvalue and position in the set that is relevant for renormalizing bit
%find indices of b0s
b0i=find([corrected_dkis.bval]<51);
% b1000i=find(abs([corrected_dkis.bval]-1000)<100);
% b2000i=find(abs([corrected_dkis.bval]-2000)<100);
% b01i=b0i(b0i<=(b1000i(end)+1));%b0s from the 1000s set
% b02i=b0i(b0i>=(b2000i(1)-1));%b0s from the 1000s set; expression is not really that general
%% coregistering (done by eddy)

% b0i(end+1)=b2000i(end)+1;
% for i=2:length(b0i)-1
%     [dkis(b0i(i):(b0i(i+1)-1)),M{i}]=coregister_obj(dkis(b0i(1)),dkis(b0i(i):(b0i(i+1)-1)));
% end
%above, i'm grabbing each b0, coreging to the first, and then applying that
%to each dwi until the next b0
%if registrations look good, i should reg everything after, too.

% get b0 average
b0s=at_read_vols(corrected_dkis(b0i)); %grab the b0s
b0sm=mean(b0s,4);
corrected_dkis(b0i(1)).img=b0sm; %first b0 replaced with mean
% later all the other b0s will be deleted



for ji=1:length(corrected_dkis)
    corrected_dkis(ji).img(isnan(corrected_dkis(ji).img)|isinf(corrected_dkis(ji).img))=0;
end %you might get errors from ratio if you don't do this

%delete other b0s 
corrected_dkis(b0i(2:end))=[]; %do this at the end so that you don't have to go back and find all the indices again

flags=[];
corrected_dkis(1).bvec=[]; %i believe this should write correctly

d2n2s_write(corrected_dkis,[basefolder filesep 'output'],'dkiss',flags)

%% RESLICE THE IMAGES 

% target_fn='/Users/andrew/re/coop/server_copy/datasets/PD/PD_MRIs2/PD_MRIs/VaP50189/dke/4D_f.nii';
% source_fn=[basefolder filesep 'output' filesep 'dkiss.nii'];
% moved_source_fn=[basefolder filesep 'output' filesep 'rdkiss.nii'];
% QC_folder='/Users/andrew/re/QC_dump';
% spm_reslice_dont_move_data(target_fn,source_fn,moved_source_fn,QC_folder)

%% make custom parameters from template

save([basefolder filesep 'workspace.mat'])
delete(gcp('nocreate'))

ndirr=[30,30];

toc
tic

dke_options('studydir',[basefolder filesep 'output'],'preprocess_options_fn_nii','dkiss.nii','ndir',ndirr,'fn_gradients',[basefolder filesep 'output' filesep 'dkiss.bvec'],'idx_gradients',{1:ndirr 1:ndirr},'T',200,'fwhm_img',1.25*[3.0 3.0 3.0])
%try dke_options('studydir',[basefolder filesep 'output' filesep 'spm_reslice'],'preprocess_options_fn_nii','redkiss.nii','ndir',ndirr,'fn_gradients',[basefolder filesep 'output' filesep 'spm_reslice' filesep 'redkiss.bvec'],'idx_gradients',{1:36 1:50},'T',400,'fwhm_img',1.25*[3.0 3.0 3.0])
%try dke_options('studydir',[basefolder filesep 'output' filesep 'shen_reslice'],'preprocess_options_fn_nii','redkiss.nii','ndir',ndirr,'fn_gradients',[basefolder filesep 'output' filesep 'shen_reslice' filesep 'redkiss.bvec'],'idx_gradients',{1:36 1:50},'T',400,'fwhm_img',1.25*[3.0 3.0 3.0])



delete(gcp('nocreate'))
study_dirr=[basefolder filesep 'output'];
dke_ft('/Users/andrew/re/coop/ft_parameters.txt',study_dirr)

toc
end













%% functions used
function [a,b]=just_eddy(dwi_dir,outdir,pick)
% use d2n2s to run eddy from 2 dwi objects
% this assumes you have dwis in a folder ALONE and topup b0s acquired differently
% in another folder ALONE

%notes: this function will not work if you have spaces in your pathname

%% handle variable input conditions
% THIS SECTION UNDER CONSTRUCTION
%make output dir if it doesn't exist
if ~exist(outdir,'dir')
    mkdir(outdir)
end

%assign to flags for reading
dflags.pick=[];

if exist('pick','var') && isfield(pick,'dwis') % could get a list of fieldnames and see if any of the strings contain 'dwi' lol
    dflags.pick=pick.dwis;
end

%% read data
% dwi_dir='';
% topup_dir='';
dwis=d2n2s(dwi_dir,dflags); %in the future, i might run flags.no='img'. This would necessitate taking dwis.bval==... and putting it into fslroi. It would also allow one to run this script from .nii.gzs

%% make a mask
%find the first b0
idx_b0s=find([dwis.bval]<51); %looking for b0s here, but some datasets are really annoying
idx_b0s_1st=idx_b0s(1);

%write to a different folder
b01=dwis(idx_b0s_1st);
intermediate_dir=[outdir filesep 'intermediate'];
mkdir(intermediate_dir)
d2n2s_write(b01,intermediate_dir,'b0',[])

%use bet
in=[intermediate_dir filesep 'b0.nii'];
mask_fn=[intermediate_dir filesep 'b0_mask.nii'];
command=['bet ' in ' ' in ' -m -n']; % these options will create a brain mask, but no skull stripped b0
[a,b]=system(command); 

%% prepare things for topup proper
%make a topup directory
topup_interm_dir=[outdir filesep 'intermediate']; %this is the same as the other intermediate dir because it's the quickest way to appropriate this script for just eddy
mkdir(topup_interm_dir)

%% use the simple_acqp function to make index and acqp text files
index_fn=[topup_interm_dir filesep 'index.txt'];
acqp_fn=[topup_interm_dir filesep 'acqp.txt'];
simple_acqp_for_just_eddy(dwis,topup_interm_dir); %this command makes a two line acqp and a corresponding index file

%% get bvec, bval, nii files for main 4d image
%pretty sure no other bvec,bval files have been written to this main folder
bvecs=dir([dwi_dir filesep '*.bvec']);
bvecs_fn=[dwi_dir filesep bvecs.name];%hope dir only grabbed one file
ensure_bvecs_are_fsl_convention(bvecs_fn)

bvals=dir([dwi_dir filesep '*.bval']);
bvals_fn=[dwi_dir filesep bvals.name]; %only one plz

dwi_image_dir=dir([dwi_dir filesep '*.nii']);
dwi_image_fn=[dwi_dir filesep dwi_image_dir.name]; %only one plz


%% run the FSL commands!
eddy_dir=[outdir filesep 'eddy'];
out_fn=[eddy_dir filesep 'eddyd'];
mkdir(eddy_dir)

command=['eddy --imain=' dwi_image_fn '  --mask=' mask_fn ' --index=' index_fn ...
    ' --acqp=' acqp_fn ' --bvecs=' bvecs_fn ' --bvals=' bvals_fn ' --out=' ...
    out_fn ' --repol --ol_nstd=20 --data_is_shelled']; % --slm=linear
[a,b]=system(command)



% function simple_acqp(json_fn,out_dir,num_imgs)
function simple_acqp_for_just_eddy(dwi_obj,out_dir)
%purpose: build acqp.txt and index.txt

%% find out how many b0s you had from each scan
%this sort of duplicates efforts from the main body of the script...
idx_b0_dwi=find([dwi_obj.bval]<51);
num_dwi_b0s=numel(idx_b0_dwi);


%% other stuff
j_struct{1}=dwi_obj(1).json; %json1);

%see if you have the fields you need
for i=1
    if ~isfield(j_struct{i},'TotalReadoutTime') || ~isfield(j_struct{i},'PhaseEncodingDirection')
        error('can''t find the required fields')
    end
end

%build two lines for acqp.txt
stringss={'j','j-'}; %known direction strings
vectorss={[0 1 0],[0 -1 0]};% corresponding vectors; 

for i=1
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
for iii=1:num_dwi_b0s %1/29/2019 AT: I think you only need the unique lines? Maybe...
    fprintf(fidd,'%i %i %i %f\n',line{1}');
end
% for jjj=1:num_top_b0s
%     fprintf(fidd,'%i %i %i %f\n',line{2}');
% end
fclose(fidd);
    


%also make index file
fidd = fopen([out_dir filesep 'index.txt'],'w');

index_array1=repmat(1,[numel(dwi_obj),1]); %this just writes as many 1s as there are dwis
fprintf(fidd,'%i\n',index_array1');

fclose(fidd);

end
end
function ensure_bvecs_are_fsl_convention(bvecs_fn)
try
    bvecs=load(bvecs_fn);
    if size(bvecs,2)==3 && size(bvecs,1)~=3 %if they're cbi convention && not ambiguous
        bvecs2=num2cell(bvecs',1); 
    elseif size(bvecs,1)==3 && size(bvecs,2)~=3 %if they're fsl/dcm2niix convention && not ambiguous
        return
    elseif size(bvecs,1)==3 && size(bvecs,2)==3 %if they're ambiguous
        return
    else 
        error('none of your bvec matrix dimensions are 3, unable to parse')
    end
    good_looking_bvecs=cat(2,bvecs2{:});
    dlmwrite(bvecs_fn,good_looking_bvecs,' ')
catch 
    warning('couldn''t ensure bvecs are fsl convention')
end
end


