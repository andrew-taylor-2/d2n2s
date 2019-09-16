function [out_fn,mask_fn,bvals_fn]=topup_n_eddy(dwi_dir,topup_dir,outdir,pick)
% use d2n2s to run eddy from 2 dwi objects
% this assumes you have dwis in a folder and topup b0s acquired differently
% in another folder

%notes: this function will not work if you have spaces in your pathname

%% handle variable input conditions

%make output dir if it doesn't exist
if ~exist(outdir,'dir')
    mkdir(outdir)
end

%% if a pick flag exists, use this to choose niftis inputs in a folder with multiple niftis
%Would love to find a way to get rid of "pick..."

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

%% make a mask
%find the first b0
idx_b0s=find([dwis.bval]<101);
idx_b0s_1st=idx_b0s(1);

%also do this for topups for later
if isfield(topups,'bval')
    idx_b0s_topups=find([topups.bval]<101);
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
simple_acqp(dwis,topups,topup_interm_dir); %this command makes a two unique line acqp and a corresponding index file

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

% get filenames from d2n2s struct
bvecs_fn=dwis(1).fns.bvec;
bvals_fn=dwis(1).fns.bval;
dwi_image_fn=dwis(1).fns.nii;

%% run the FSL commands!
eddy_dir=[outdir filesep 'eddy'];
out_fn=[eddy_dir filesep 'eddyd'];
mkdir(eddy_dir)


command=['topup --imain=' topup_b0s_fn ' --datain=' acqp_fn ' ' ...
    '--config=b02b0.cnf --out=' topup_interm_dir filesep 'my_topup_results ' ...
    '--fout=' topup_interm_dir filesep 'my_field --iout=' topup_interm_dir filesep ...
    'my_unwarped_images'];
disp(command)
system(command)

command=['eddy --imain=' dwi_image_fn '  --mask=' mask_fn ' --index=' index_fn ...
    ' --acqp=' acqp_fn ' --bvecs=' bvecs_fn ' --bvals=' bvals_fn ' --out=' ...
    out_fn ' --topup=' topup_interm_dir filesep 'my_topup_results --repol --data_is_shelled'];
disp(command)
system(command)



% function simple_acqp(json_fn,out_dir,num_imgs)
function simple_acqp(dwi_obj,topup_obj,out_dir)
%purpose: build acqp.txt and index.txt

%% find out how many b0s you had from each scan
%this sort of duplicates efforts from the main body of the script...
idx_b0_dwi=find([dwi_obj.bval]<101);
num_dwi_b0s=numel(idx_b0_dwi);

if isfield(topup_obj,'bval')
    idx_b0_topups=find([topup_obj.bval]<101);
else %assume they're all b0s
    warning('couldn''t find bvals for topup images -- assuming they''re all b0s')
    idx_b0_topups=1:numel(topup_obj);
end
num_top_b0s=numel(idx_b0_topups);

%% parse json info, match to known PE dir strings
% j_struct{1}=readJson(json_fn); 
j_struct{1}=dwi_obj(1).json; 
j_struct{2}=topup_obj(1).json;

%see if you have the fields you need
for i=1:2
    if ~isfield(j_struct{i},'TotalReadoutTime') || ~isfield(j_struct{i},'PhaseEncodingDirection')
        error('can''t find the required fields')
    end
end

%build two lines for acqp.txt
stringss={'j','j-','i','i-'}; %known direction strings
vectorss={[0 1 0],[0 -1 0],[1 0 0],[-1 0 0]};% corresponding vectors;

for i=1:2
    time=j_struct{i}.TotalReadoutTime;
    dir_string{i}=j_struct{i}.PhaseEncodingDirection;
    
    match_bool=strcmp(stringss,dir_string{i});
    assert(sum(match_bool)==1,'Couldn''t read your PEdirection string')
    
    %match vector representation
    vec=vectorss{match_bool};
    
    %this is a unique line of acqp.txt
    line{i}=[vec time]; 
end


%% write .txt files 
% there's one unique line for dwis and one for topups; this is going to
% write them as many times as each set has b0s (bc that's what topup wants)

% andrew sept, 2019: am I sure we need the unique lines to be repeated?
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

index_array1=repmat(1,[numel(dwi_obj),1]); %this just writes as many 1s as there are dwis -- the one tells eddy that the first line of acqp.txt is how your --imain images were acquired
fprintf(fidd,'%i\n',index_array1');

%you might use the next 3 lines if you were to append the topup b0s to the
%end of the dwi set
% index_of_second_unique_acqp_line=num_top_b0s+1;
% index_array2=repmat(index_of_second_unique_acqp_line,[numel(topup_obj),1]); %this just writes as many of the index of the second unique acqp line as there are topups
% fprintf(fidd,'%i\n',index_array2');

fclose(fidd);

end

end

function ensure_bvecs_are_fsl_convention(bvecs_fn)
try
    bvecs=load(bvecs_fn);
    if size(bvecs,2)==3 && size(bvecs,1)~=3 %if they're cbi convention && not ambiguous
        bvecs2=num2cell(bvecs',1); 
    elseif size(bvecs,1)==3 && size(bvecs,2)~=3 %if they're fsl/dcm2niix convention && not ambiguous
%         bvecs2=num2cell(bvecs,1);
        return
    elseif size(bvecs,1)==3 && size(bvecs,2)==3 %if they're ambiguous
%         bvecs2=num2cell(bvecs,1);
%         warning('your bvecs are square, so the orientation convention is ambiguous to this function. assuming they are in fsl/dcm2niix convention. user may need to transpose.')
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

%% begin external code from Chris Rorden's nii_preprocess
%this is unused as of now, but sometimes half sphere bvec acquisitions
%benefit greatly from using the --slm=linear option in eddy

function isHalfSphere = HalfSphere(bvec_nam)
%return true if DWI vectors sample half sphere
% bvec_nam : filename of b-vector file
%Example
% HalfSphere('dki.bvec')
% if ~exist('bvec_nam','var') %file not specified
%    fprintf('Select b-vec file\n');
%    [A,Apth] = uigetfile({'*.bvec'},'Select b-vec file');
%    if isnumeric(A), return; end; %nothing selected
%    bvec_nam = strcat(Apth,char(A));
% end;
% %handle different extensions bvec/nii/nii.gz
[p,n,x]=fileparts(bvec_nam);
if ~strcmpi(x,'.bvec')
    fnm = fullfile(p,[n,'.bvec']);
    if ~exist(fnm,'file') %assume img.nii.gz
        [~,n,x] = fileparts(n);
        fnm = fullfile(p,[n,'.bvec']);
    end
    if ~exist(fnm,'file')
        error('Unable to find bvec file for %s', vNam);
    end
    bvec_nam = fnm;
end
isHalfSphere = false;
bvecs = load(bvec_nam); 
bvecs = unitLengthSub(bvecs)';
bvecs(~any(~isnan(bvecs')),:) = [];
if isempty(bvecs), return; end;
mn = unitLengthSub(mean(bvecs)')';
if isnan(mn), return; end;
mn = repmat(mn,size(bvecs,1),1);
minV = min(dot(bvecs',mn'));
thetaInDegrees = acosd(minV);
if thetaInDegrees < 110
    isHalfSphere = true;
    fprintf('Sampling appears to be half sphere (%g-degrees): %s\n', thetaInDegrees, bvec_nam);
end;
end

function Xu = unitLengthSub(X) 
%Use pythagorean theorem to set all vectors to have length 1
%does not require statistical toolbox 
if size(X,1) ~= 3, error('expected each column to have three rows (X,Y,Z coordinates)'); end;
Dx =  sqrt(sum((X.^2)));
Dx = [Dx; Dx; Dx];
Xu = X./Dx;
end 
%% end external code