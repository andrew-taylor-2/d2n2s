%% RDF paper processing 
clear all
clc

%% to do 
% make 1 dicom folder with all dki images

%% Convert DKI to nifti 

dke_exe='C:\programs\DKE_2.6.0i\Executables\dke ';
dke_ft_exe='C:\programs\DKE_Tracktography\dke_ft.exe ';
parentFolder='V:\taylor_j\test_folder\6283\';
subj_list={'6283'};
    
for i=1:length(subj_list)
    b12root = [parentFolder subj_list{i}];
    list=dir(b12root);
    list4=dir([list(1).folder '\dicom']);
    %make 4D niis from main dwi dicoms
    dwihdr=dicominfo([list(1).folder '\dicom' '\' list4(3).name]);
    grad=dcm2nii_dki(list(1).folder,lower(dwihdr.SeriesDescription));

%     seqroot = [b12root '/nifti/' seq{k}]
    niftiNames=dir([b12root '\nifti\*.nii']);

%     
% %% TOPUP
% %  make sure this is in the right place (ie there is no interpolation or
% %  whatever that would mess with other functions you use
%     
%     %locate topup dicoms
%     d=dir([b12root '\*TOPUP*']);
%     if length(d)==1
%         topupdir=[d(1).folder '\' d(1).name];
%     elseif length(d)==2
%         topupdir=[d(2).folder '\' d(2).name];
%     end
% %    topup_dicomsPA='/dicom/FBWM_DKI_b0_PA_TOPUP_9/dicom';
% 
%     %get topup-relevant info from dicoms
%     d=dir([topupdir]);
%     dcm_hdr=dicominfo([topupdir '\' d(3).name]);
%     Effective_echo_spacing_PA=1/(dcm_hdr.Private_0019_1028*double(dcm_hdr.AcquisitionMatrix(4)));
%     Total_read_out_PA=Effective_echo_spacing_PA*(dcm_hdr.NumberOfPhaseEncodingSteps-str2num(dcm_hdr.Private_0051_1011(5)));
% 
%     Effective_echo_spacing_AP=1/(dwihdr.Private_0019_1028*double(dwihdr.AcquisitionMatrix(4)));
%     Total_read_out_AP=Effective_echo_spacing_PA*(dwihdr.NumberOfPhaseEncodingSteps-str2num(dwihdr.Private_0051_1011(5)));
% 
%     %make nii from topup dicoms
%     mkdir([topupdir '\dicom']);
%     source=[topupdir '\*.dcm'];
%     target=[topupdir '\dicom'];
%     movefile(source,target);
% 
%     mkdir([topupdir '\nifti']);
%     copyfile([b12root '\nifti\' niftiNames(1).name], [topupdir '\nifti']);
%     movefile([topupdir '\nifti\' niftiNames(1).name],[topupdir '\nifti\' 'temp.nii']);
%     grad=dcm2nii_dki_b0_vitria(topupdir,lower(dcm_hdr.SeriesDescription));
%      
%     mkdir([b12root '\b0s']);
%     copyfile([topupdir '\nifti\img000_b0_01.nii' ],[b12root '\b0s\' 'topb0.nii']);
%     copyfile([b12root '\nifti\img001_ep_b0.nii'], [b12root '\b0s\' 'firstb0.nii']);
%     listb0s=dir([b12root '\nifti\img*_ep_b0.nii']);
%     for i1=2:length(listb0s)
%         copyfile([b12root '\nifti\' listb0s(i1).name], [b12root '\b0s\' 'endb0_' listb0s(i1).name(6) '.nii']);
%     end
%     %all b0s should now be in a folder called 'b0s'
%     
%     %make all_b0s file
%     olist=dir([b12root '\b0s\end*']);
%     make_4D_nii([b12root '\b0s'], {'firstb0.nii' olist(1:end).name 'topb0.nii'},'all_b0.nii');
%     
% 
%     %dicm2nii scripts will be needed in this section (and, if I revise,
%     %possibly other sections
% 
%     
% %     files_AP=dir([]);
% %     files_PA=dir([12root '/gibbsringing/' 'topup*b0*.nii']);
% % 
% %     make_4D_nii([b12root '/gibbsringing/'],{files_AP.name files_PA.name},'4D_all_b0.nii');  
% %     mkdir([b12root '/topup']);
% %     movefile([b12root '/gibbsringing/4D_all_b0.nii'],[b12root '/topup/4D_all_b0.nii'])
%     text_file_topup=[0 -1 0 Total_read_out_AP];
%     text_file_topup=repmat(text_file_topup,[length(olist)+1,1]);
%     text_file_topup=[text_file_topup; repmat([0 1 0 Total_read_out_PA],[1,1])]; % in that last vector, the first element, "1", refers to the fact that we have one topup image. maybe this will not be assumed in later versions of this script
%     save([b12root '/b0s/acq_parameters.txt'],'text_file_topup','-ascii');
%     
% 
%     command=['topup --imain=' b12root '/b0s/all_b0.nii --datain=' b12root '/b0s/acq_parameters.txt ' '--config=b02b0.cnf ' ' --out=' b12root '/b0s/topup_results ' ' --fout=' b12root '/b0s/the_field ' ' --iout=' b12root '/b0s/corrected_b0s '];
%     [~,cmdout] = system(command);
%     
%     %examine "corrected b0s" to make sure they look good and/or better than
%     %the b0s did before
%     
%     make_4D_nii([b12root '\nifti\'],{niftiNames.name},'4D_AP.nii');
%     
%     command=['applytopup --imain=' b12root '/nifti/4D_AP.nii --datain=' b12root '/b0s/acq_parameters.txt ' '--inindex=1 ' ' --topup=' b12root '/b0s/topup_results ' ' --out=' b12root '/b0s/4Dc ' ' --method=jac'];
%     [~,cmdout] = system(command);    
%     

    
    
%% Denoise
%           CHANGE THESE TO THE ONES ON YOUR MACHINE    
command = sprintf('C:/msys64/mingw64 C:/msys64/home/mrtrix3/release/bin/dwidenoise %s/4Dc.nii %s/4D_DN.nii -noise %s/noise.nii', [b12root '/b0s'], [b12root '\nifti'], [b12root '\nifti']);
[status,cmdout] = system(command);

disp('MATLAB continues after calling EXE')
% Check to see if the EXE process exists
flag = true;
while flag
     disp('EXE still running');   
     flag = isprocess('mintty.exe');
     pause(1) 
end
disp('EXE Done')
disp('continuing with MATLAB script')

command = sprintf('C:/msys64/mingw64 C:/msys64/home/mrtrix3/release/bin/mrcalc %s/4Dc.nii %s/4D_DN.nii -subtract %s/res.nii', [b12root '\b0s'], [b12root '\nifti'], [b12root '\nifti']);
[status,cmdout] = system(command);

disp('MATLAB continues after calling EXE')
% Check to see if the EXE process exists
flag = true;
while flag
     disp('EXE still running');   
     flag = isprocess('mintty.exe');
     pause(1) 
end
disp('EXE Done')
disp('continuing with MATLAB script')


%% unring dki images 
    DN=spm_read_vols(spm_vol(fullfile([b12root '\nifti'],'4D_DN.nii')));
    list=dir(fullfile([b12root '\nifti'],'img*'));
    [dim1,dim2,dim3,dim4]=size(DN);
    parfor j=1:dim4
    img(:,:,:,j)=unring(DN(:,:,:,j));
    hdr=spm_vol(fullfile([b12root '\nifti'],[list(j).name]));
    hdr.dt=[16 0];
    int=img(:,:,:,j);
    int(isnan(int))=0;
    spm_write_vol(hdr,int);
    end

    end
    
 %% Coregister AVE2 to AVE1 b0
fprintf('Co-registering images...\n')
fn_target = fullfile(b12root, '/nifti/ave1/img001_ep_b0.nii');
fn_source = fullfile(b12root, '/nifti/ave2/img001_ep_b0.nii');   % source file is the first b = 0 image in the series returned by the operating system
M=coregister(fn_target, fn_source, fullfile(b12root, 'nifti/ave2'),'.nii');

fprintf('Co-registration complete.\n')

mkdir(fullfile(b12root, 'nifti/ave2_coreg'));
list=dir(fullfile(b12root, 'nifti/ave2/rimg*'));
for j=1:length(list)
    fn_source=fullfile(b12root, 'nifti/ave2',[ list(j).name ]);
    fn_target=fullfile(b12root, 'nifti/ave2_coreg',[ list(j).name ]);
    movefile(fn_source,fn_target)
end

% Coregister b0 to AVE1 b0
fprintf('Co-registering images...\n')
fn_target = fullfile(b12root, '/nifti/ave1/img001_ep_b0.nii');
fn_source = fullfile(b12root, '/nifti/b0/img000_b0_01.nii');   % source file is the first b = 0 image in the series returned by the operating system
M=coregister(fn_target, fn_source, fullfile(b12root, 'nifti/b0'),'.nii');

fprintf('Co-registration complete.\n')

mkdir(fullfile(b12root, 'nifti/b0_coreg'));
list=dir(fullfile(b12root, 'nifti/b0/rimg*'));
for j=1:length(list)
    fn_source=fullfile(b12root, 'nifti/b0',[ list(j).name ]);
    fn_target=fullfile(b12root, 'nifti/b0_coreg',[ list(j).name ]);
    movefile(fn_source,fn_target)
end    

%% average b0's

mkdir(fullfile(b12root, 'nifti/combined'));
list = dir(fullfile(b12root, 'nifti/b0_coreg', '*_b0_*.nii'));
hdr = spm_vol(fullfile(b12root, 'nifti/b0_coreg', list(1).name));
imgavg = spm_read_vols(hdr);
for j = 2:length(list)
    hdr = spm_vol(fullfile(b12root,'nifti/b0_coreg',list(j).name));
    img = spm_read_vols(hdr);
    imgavg = imgavg + img;
end
hdr = spm_vol(fullfile(b12root, 'nifti/ave1','img001_ep_b0.nii'));
img = spm_read_vols(hdr);
imgavg = imgavg + img;
hdr = spm_vol(fullfile(b12root, 'nifti/ave2_coreg','rimg001_ep_b0.nii'));
img = spm_read_vols(hdr);
imgavg = imgavg + img;
imgavg = imgavg / (length(list)+2);

hdr.dt=[16 0];
hdr.fname = fullfile(b12root, 'nifti/combined', 'b0_avg.nii');
imgavg(isnan(imgavg))=0;
spm_write_vol(hdr, imgavg);

end

%% AVERAGE  AVE1 &2

list = dir(fullfile(b12root, 'nifti/ave1', '*b1000.nii'));

for j = 1:length(list)
    hdr = spm_vol(fullfile(b12root, 'nifti/ave1', list(j).name));
    imgavg = spm_read_vols(hdr);
    hdr = spm_vol(fullfile(b12root,'nifti/ave2_coreg', ['r' list(j).name]));
    img = spm_read_vols(hdr);
    imgavg = imgavg + img;
    imgavg = imgavg / (2);
    hdr.dt=[16 0];
    hdr.fname = fullfile(b12root, 'nifti/combined', list(j).name);
    imgavg(isnan(imgavg))=0;
    spm_write_vol(hdr, imgavg);
end

list = dir(fullfile(b12root, 'nifti/ave1', '*b2000.nii'));

for j = 1:length(list)
    hdr = spm_vol(fullfile(b12root, 'nifti/ave1', list(j).name));
    imgavg = spm_read_vols(hdr);
    hdr = spm_vol(fullfile(b12root,'nifti/ave2_coreg', ['r' list(j).name]));
    img = spm_read_vols(hdr);
    imgavg = imgavg + img;
    imgavg = imgavg / (2);
    hdr.dt=[16 0];
    hdr.fname = fullfile(b12root, 'nifti/combined', list(j).name);
    imgavg(isnan(imgavg))=0;
    spm_write_vol(hdr, imgavg);
end

%% make 4D nifti and run DKE and FT
%make list of img* images and then strip out the ones from the end, using
%the length of olist
files = dir([b12root '\nifti\' 'img*.nii']);
make_4D_nii([b12root '\nifti\'],{files(1:(length(files)-length(olist))).name},'4d.nii');
mkdir([b12root '\dke\'])
movefile([b12root '\nifti\4d.nii'],[b12root '\dke\4d.nii'])
 img=spm_read_vols(spm_vol([b12root '\dke\4d.nii']));
 img(isnan(img))=0;
 make_4D_nii(spm_vol([b12root '\dke\4d.nii']),img,'4d.nii');

mkdir([b12root '/output']);
movefile([b12root '\dke_parameters.txt'],[b12root '\dke'])
command=[dke_exe b12root '/dke/dke_parameters.txt '];
[status,cmdout] = system(command,'-echo');

movefile([b12root '\ft_parameters.txt'],[b12root '\dke'])
command=[dke_ft_exe b12root '/dke/ft_parameters.txt '];
[status,cmdout] = system(command,'-echo');
