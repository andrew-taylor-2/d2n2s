%common paths
dke_path = 'C:\Programs\DKE\win64_internal\Executables\dke';        %full path to dke.exe file
ft_path = 'C:\Programs\DKE_Tracktography\dke_ft.exe';       %full path to dke_ft.exe file
mingw64_path = 'C:\Programs\msys64\mingw64.exe';                         %path to mrtrix3 mingw terminal (blue)
denoise_path = 'C:\Programs\msys64\home\Andrew\mrtrix3\bin';
dcm2nii_path='C:\Programs\MRIcroGL\dcm2niix.exe';

%prior to this script, dicoms should be sorted and converted with dcm2niix

% dicom_sort('C:\re\test\AD_study','')
% eval(sprintf('!dcm2niix %s',folder))

fbi_dir='C:\re\test\AD_study\dicom\FBI_b6000_128_16';%spm_select(Inf,'dir','choose fbi dir');
dki_dir='C:\re\test\AD_study\dicom\DKI_30_Dir_17';%spm_select(Inf,'dir','choose dki dir');


command=[dcm2nii_path ' ' fbi_dir];
[a,b]=system(command)
command=[dcm2nii_path ' ' dki_dir];
[a,b]=system(command)




%d2n2s isn't well suited to putting things into mrtrix, so
%i'll just denoise first and then read in my images.

%denoise
fbinii=dir([fbi_dir '/*.nii']);
[status,cmdout]=mrtrix_denoise([fbinii.folder filesep fbinii.name],'_DN',mingw64_path,denoise_path);

dkinii=dir([dki_dir '/*.nii']);
[status,cmdout]=mrtrix_denoise([dkinii.folder filesep dkinii.name],'_DN',mingw64_path,denoise_path);


dkdn=dir([dki_dir '/*_DN.nii']);
flags1.pick=[dkdn.folder filesep dkdn.name];
dkis=d2n2s(dki_dir,flags1);

fbdn=dir([fbi_dir '/*_DN.nii']);
flags2.pick=[fbdn.folder filesep fbdn.name];
fbis=d2n2s(fbi_dir,flags2);



%store one from each dataset to compare unringing results
examples=join_obj(dkis(1),fbis(1));

dkis=unringfun(dkis);
fbis=unringfun(fbis);


%coreg dkis to fbi
[rdkis,M1]=coregister_obj(fbis(1),dkis);
% M1=estimate_reg_obj(fbis(1),dkis);

for i=1:length(rdkis)
    rdkis(i).bvec=M1(1:3,1:3)*rdkis(i).bvec;
end

%or this
% a=[rdkis.bvec];
% a=M1(1:3,1:3)*a;
% a=num2cell(a,1);
% [rdkis.bvec]=a{:};


%average b0s:
%find indices of b0s and "read" images
%in all previous tests, using at_read_vols() is the same as just
%using the .img data, but we're doing this just to be safe
b0i=find([rdkis.bval]==0);
b0s=cat(4,rdkis(b0i).img);

b0i2=find([fbis.bval]==0);
b0s2=cat(4,fbis(b0i2));

allb0s=cat(4,b0s,b0s2);

b0sm=mean(allb0s,4);

rdkis(b0i(1)).img=b0sm; %first b0 replaced with mean
rdkis(b0i(2:end))=[]; %other b0s deleted
fbis(b0i2)=[]; %other b0s deleted


%change dke_parameters and write stuff
mkdir([dki_dir '\dke'])
mkdir([fbi_dir '\final'])
flags=[];
d2n2s_write(fbis,[fbi_dir '\final'],'fbis',flags)
d2n2s_write(rdkis,[dki_dir '\dke'],'dkis',flags)

command = ['C:\Programs\DKE\win64_internal\Executables\dke.exe ' dki_dir '\dke\dke_params.txt'];
[status,cmdout] = system(command);

optimize_FBI([fbi_dir '\final'])