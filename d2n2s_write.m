function d2n2s_write(dwi,folder,name,flags)
%take matlab dwi object and write it to file
% input is "object" (structure created by d2n2s), output is files
%from andrew's OOP scripts
%flags
%  flags.gr - 1 for writing the cbi way. default is 1. 0 for the fsl way
%  flags.vol - 3 for 3d files. 4 for a 4d file
%  flags.nfn - 1 for writing 3d volumes according to 'folder','name'.
%            - 0 for using the 'dwi.fn' field to write. 1 is default

%if you wanted to use this on only one field, you could make another object
%with just the field. or use a segment of the object

if ~exist(folder,'dir')
    try
        mkdir(folder)
    catch
    end
end

if isfield(dwi,'bval')
    bvalw=[dwi(:).bval];
    fn=fullfile(folder,[name '.bval']);
    fid=fopen(fn,'w');
    fprintf(fid,'%u ',bvalw);
    fclose(fid);
end

if isfield(dwi,'bvec')
    bvecs=[dwi(:).bvec];
    fn=fullfile(folder,[name '.bvec']);
    if ~isfield(flags,'gr') || isempty(flags.gr)
        flags.gr=1;
    end
    if flags.gr
        if size(bvecs,1)>size(bvecs,2) %they're cbi convention
        elseif size(bvecs,1)<size(bvecs,2) %they're fsl convention
            bvecs=bvecs';
        end
        dlmwrite(fn,bvecs,' ');
    elseif ~flags.gr
        if size(bvecs,1)>size(bvecs,2) %they're cbi convention
            bvecs=bvecs';
        elseif size(bvecs,1)<size(bvecs,2) %they're fsl convention
        end
        dlmwrite(fn,bvecs,' ')
    end
end

if isfield(dwi,'hdr') && isfield(dwi,'img')
    fn=fullfile(folder,[name '.nii']);
    if ~isfield(flags,'vol') || isempty(flags.vol) || ~any(flags.vol==[3,4])
        flags.vol=4;
    end
    if flags.vol==4
        fn=fullfile(folder,[name '.nii']); %this is named according to the function's inputs
        if exist(fn,'file'); warning('you''re writing to an image file that already exists; note that if you overwrite an existing image with an image containing less volumes, the extra volumes at the end will stay');end
        for i1=1:length(dwi)
            dwi(i1).hdr.fname=fn;
            dwi(i1).hdr.n=[i1,1];
            spm_write_vol(dwi(i1).hdr,dwi(i1).img);
        end
    end
    if flags.vol==3
        if ~isfield(flags,'nfn') || isempty(flags.nfn) || ~any(flags.vol==[0,1])
            flags.nfn=1;
        end
        if flags.nfn==0
            if ~isfield(dwi,'fn')
                error('you''re going to need new filenames to write stuff.')
            end
            for i2=1:length(dwi)
                dwi(i2).hdr.fname=dwi(i2).fn; %the volumes are named according to the dwi object .fn field
                spm_write_vol(dwi(i2).hdr,dwi(i2).img)
            end
        end
        
        if flags.nfn==1
            for i2=1:length(dwi)
                fn=fullfile(folder,[name sprintf('_%05d.nii',i2)]); %the volumes are programmatically named, similar to how
                %spm_split_file would name them
                dwi(i2).hdr.fname=fn;
                dwi(i2).hdr.n=[1,1];
                spm_write_vol(dwi(i2).hdr,dwi(i2).img)
            end
        end
    end
end
if isfield(dwi,'json') 
    try; if numel(dwi)>1 && ~isequal(dwi.json)
        warning('Not every json in the dataset to be written is equal to one another. Writing only the first element''s .json field')
    end; catch; end
    jsontext=jsonencode(dwi(1).json);
    fn=fullfile(folder,[name '.json']);
    fid=fopen(fn,'w');
    fprintf(fid,'%s',jsontext);
    fclose(fid);
end
