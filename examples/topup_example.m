%% new example script
%% new script for lorna based off:
%preprocessing with FSL only
% -andrew taylor, 10/16

%% use dcm2niix, make sure it outputs .json files

%prior to this script, dicoms should be sorted and converted with dcm2niix
%e.g.:
% dicom_sort('C:\re\test\AD_study','')
% eval(sprintf('!dcm2niix %s',folder))

% to be edited by the user
dki_dir='/Users/andrew/re/test/topup_n_eddy/20180323_084032_Test_ADStudy/dicom/DKI_30_Dir_17';
% fbi_dir='/Users/andrew/re/test/topup_n_eddy/20180323_084032_Test_ADStudy/dicom/FBI_b6000_128_16';
top_dir='/Users/andrew/re/test/topup_n_eddy/20180323_084032_Test_ADStudy/dicom/DKI_30_Dir_TOP_UP_PA_23';
out_dir='/Users/andrew/re/test/topup_n_eddy/20180323_084032_Test_ADStudy/dicom/out';

%Notes: 
%The dki_dir should contain a .nii, .bvec, .bval, and .json file.
%The top_dir should contain a .nii, .json, and optionally .bvec and .bval
%files. If not, the script will assume the bvec/bval files of a b0 and just
%give you warnings you can ignore. These directories can contain any other
%file types you would like, such as dicoms, but SHOULD NOT contain other
%.nii,.bvec,.bval, or .json files. The "out_dir" doesn't have to exist yet.


if strcmp(out_dir(end),'/'); out_dir=out_dir(1:end-1);end % this is just making sure out_dir doesn't have a trailing slash since that's assumed


dkis=d2n2s(dki_dir,[]);
% fbis=d2n2s(fbi_dir,[]);

%concatenate sequences and write
% dwis=join_obj(dkis,fbis);
dwis=dkis; % if you had FBIs, you would not use this line, but rather delete it and uncomment the previous 3 lines that are actual code relating to FBI

working_dir=[out_dir filesep 'input_dwis'];
d2n2s_write(dwis,working_dir,'img',[])


%% initial corrections
dkinii=dir([working_dir '/*.nii']); %make sure it's only chosen one
fn=[dkinii.folder filesep dkinii.name];

% put fn through corrections. this function is at the bottom of this
% script.
out=DN_UR_RC(fn);

%% handle output of corrections

% try;gunzip([out '.gz']);catch;end % depends on your FSL config, but it might be spitting out .nii.gz files instead of .nii files. This should handle either case.
% actually, having both .nii and .nii.gz files will cause fsl to error.
pick.dwis=out;

%% run topup_and_eddy

% topup_n_eddy(dwi_dir,topup_dir,outdir,pick)
[final_out_fn,final_mask_fn,final_bval_fn]=topup_n_eddy(working_dir,top_dir,out_dir,pick); % the outputs are just names of files created that you'll need
%out bvecs name
eddy_out_bvecs=[final_out_fn '.eddy_rotated_bvecs'];

%% handle outputs

%get output bvecs and take transpose
cbvecs=importdata(eddy_out_bvecs);
final_bvecs=[fileparts(final_out_fn) filesep 'final_bvecs.txt']; % I'm putting the final bvecs in the output eddy folder. Seems kind of arbitrary, feel free to change it
fprintf(fopen(final_bvecs,'w'),'%f %f %f \n',cbvecs);
fclose('all');


%gunzip data
gunzip([final_out_fn '.nii.gz'])
gunzip([final_mask_fn '.gz'])
mkdir([out_dir filesep 'dke'])

%paths for DKE
path_data=[final_out_fn '.nii'];
path_gradient=final_bvecs;
path_brain_mask=final_mask_fn;
path_bval=final_bval_fn;
path_output=[out_dir filesep 'dke/']; %because of how dke_proc works, this has to have a slash at the end. Emilie wrote that one so I'm not gonna try to change it
path_dke_parameters='~/bin/bin2/FBWM/DKEParameters.txt'; %wherever her example one is

dke_proc(path_data,path_gradient,path_brain_mask,path_bval,path_output,path_dke_parameters) %dke_proc not included


function [out,a,b]=DN_UR_RC(fn)

[path,name,ext]=fileparts(fn);

a={};
b={};
%denoise
in=fn;
out=[path filesep name '_dn' ext];
noise_out=[path filesep 'noise' ext];
command=['dwidenoise ' in ' ' out ' -noise ' noise_out ' -force'];
[a{end+1},b{end+1}]=system(command);

noise_residuals_out=[path filesep 'res' ext];
command=['mrcalc ' in ' ' out ' -subtract ' noise_residuals_out];
[a{end+1},b{end+1}]=system(command);

%gibbs unring
in=[path filesep name '_dn' ext]; % or in=out;
out=[path filesep name '_dn_ur' ext];
command=['mrdegibbs ' in ' ' out ' -force'];
[a{end+1},b{end+1}]=system(command);


%rician
in=out;
out=[path filesep name '_dn_ur_rc' ext];
signal_sqr=[path filesep name '_sqr' ext];
noise_sqr=[path filesep 'noise_sqr' ext];
diff=[path filesep 'diff' ext];

c1=['fslmaths ' in ' -sqr ' signal_sqr ';'];
c2=['fslmaths ' path filesep 'noise -sqr ' noise_sqr ';'];
c3=['fslmaths '  signal_sqr ' -sub ' noise_sqr ' ' diff ';'];
c4=['fslmaths ' diff ' -sqrt ' out];
[a{end+1},b{end+1}]=system([c1 c2 c3 c4]);
end

