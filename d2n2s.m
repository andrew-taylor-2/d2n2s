function dwi=d2n2s(dcm2niixd_folder,flags)
%radical change loll
% operates on files/folders
%
%flags.pick - if this is specified, this nii filename will be chosen to have it's
%header/image information loaded. just in case, for instance, one has
%denoised after running dcm2niix

%flags.no  - tells the program not to read certain things. i'll add more
%about this later

%it also might be useful to make the dir calls a wrapper, so that d2n2s can
%be made more general

if isstring(dcm2niixd_folder)
    dcm2niixd_folder=char(dcm2niixd_folder);
end
if isfield(flags,'pick') && exist(flags.pick) && isstring(flags.pick)
    flags.pick=char(flags.pick);
end

if ~isfield(flags,'no') || isempty(flags.no) 
    flags.no='';
end

%find files
dbvec=dir([dcm2niixd_folder filesep '*.bvec']);
dbval=dir([dcm2niixd_folder filesep '*.bval']);
%dnii=dir([dcm2niixd_folder filesep '*.nii']);
%dnii gets set later and is part of a conditional
djson=dir([dcm2niixd_folder filesep '*.json']);

%warnings
if length(dbvec)>1
    warning('found more than one bvec file. using the ''first''.')
    dbvec=dbvec(1);
end

if length(dbval)>1
    warning('found more than one bval file. using the ''first''.')
    dbval=dbval(1);
end
% if length(dnii)>1
%     warning('found more than one nii file. using the ''first''.')
% end
if length(djson)>1
    warning('found more than one json file. using the ''first''.')
    djson=djson(1);
end

%assign
if ~isempty(dbval) && ~contains(flags.no,'bval','IgnoreCase',1)
    bvals=load([dbval.folder filesep dbval.name]);
    bvals2=num2cell(bvals);
    [dwi(1:length(bvals2)).bval]=bvals2{:}; %could've assigned with the same form as bvecs below, but the length of dwi must be initialized for these assignment to work
end

if ~isempty(dbvec) && ~contains(flags.no,'bvec','IgnoreCase',1)
    bvecs=load([dbvec.folder filesep dbvec.name]);
    if size(bvecs,2)==3 && size(bvecs,1)~=3 %if they're cbi convention && not ambiguous
        bvecs2=num2cell(bvecs',1); 
    elseif size(bvecs,1)==3 && size(bvecs,2)~=3 %if they're fsl/dcm2niix convention && not ambiguous
        bvecs2=num2cell(bvecs,1);
    elseif size(bvecs,1)==3 && size(bvecs,2)==3 %if they're ambiguous
        bvecs2=num2cell(bvecs,1);
        warning('your bvecs are square, so the orientation convention is ambiguous to this function. assuming they are in fsl/dcm2niix convention. user may need to transpose.')
    else 
        error('none of your bvec matrix dimensions are 3, unable to parse')
    end
    %this doesn't generally work for 3x3 bvecs
    [dwi(1:length(bvecs2)).bvec]=bvecs2{:};     
end

if ~isfield(flags,'pick') || isempty(flags.pick) || isempty(dir(flags.pick)) %if the user has not opted to specify the nii fn himself
    dnii=dir([dcm2niixd_folder filesep '*.nii']);
    
    %sanitize input:remove 'res.nii' and 'noise.nii'
    %these are common ancillary image names in our pipeline
    problem_index=[];
    problem_index2=[];
    for j=1:length(dnii)
        if strcmpi(dnii(j).name,'res.nii')
            problem_index=[problem_index j];
        end
        if strcmpi(dnii(j).name,'noise.nii')
            problem_index2=[problem_index2 j];
        end
    end
    problem_indices=union(problem_index,problem_index2);
    dnii(problem_indices)=[];
    
elseif isfield(flags,'pick') && ~isempty(flags.pick) && ~isempty(dir(flags.pick))
    dnii=dir(flags.pick);
end

if ~isempty(dnii) && ~contains(flags.no,'nii','IgnoreCase',1) && ~contains(flags.no,'img','IgnoreCase',1)
    hdr=spm_vol([dnii.folder filesep dnii.name]);
    if exist('dwi','var')
        if length(hdr)~=length(dwi)
            warning('spm_vol thinks your 4D volume has a different number of images from your bvals file. this function might not behave as intended')
        end
    end
    if length(dnii)>1
        warning('found more than one nii file. using the ''first''.')
        dnii=dnii(1);
    end
    
    
    for i=1:length(hdr)
        hdr2{i}=hdr(i);
    end
    [dwi(1:length(hdr2)).hdr]=hdr2{:};
    
    img=spm_read_vols(hdr);
    img2=num2cell(img,[1 2 3]);
    [dwi.img]=img2{:};
    
    v0=spm_file_split2([dnii.folder filesep dnii.name]);
    for i=1:size(v0,1)
        dwi(i).fn=v0(i,:);
        dwi(i).hdr.private.dat.fname=v0(i,:); %this line is included to protect the user, as assigning to dat(:) at all in matlab will overwrite whatever is contained in hdr.private.dat.fname. this is a property of the @file_array object. and it's spooky.
    end
end
if ~isempty(djson) && ~contains(flags.no,'json','IgnoreCase',1)
    jsfn=[djson.folder filesep djson.name];
    jsvals=readJson(jsfn);
    [dwi.json]=deal(jsvals);
end

%functions used

    function Vo = spm_file_split2(V, odir)  % I'm using parts of spm_file_split to output generic new filenames for convenience. But I've taken out the part where it actually writes new images to disk.
        % Convert a 4D volume file into a series of 3D volume files
        % FUNCTION Vo = spm_file_split(V, odir)
        % V           - filename or spm_vol struct
        % odir        - output directory [default: same as input]
        %
        % Vo          - spm_vol struct array of output files
        %__________________________________________________________________________
        % Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging
        
        % Guillaume Flandin
        % $Id: spm_file_split.m 4854 2012-08-22 13:29:10Z ged $
        
        if ~nargin
            [V, sts] = spm_select(1,'nifti','Select a 4D volume file to split');
            if ~sts, return; end
        end
        if ischar(V)
            [p,n,e] = spm_fileparts(V);
            V = spm_vol(fullfile(p,[n e]));
        end
        
        [p,n,e] = spm_fileparts(V(1).fname);
        if nargin < 2
            if isempty(p), p = pwd; end
            odir = p;
        end
        
        Voo = cell(numel(V),1);
        %spm_progress_bar('Init',numel(V),'Splitting 4D Volume','Volumes Complete');
        for i1=1:numel(V)
            Voo{i1}        = fullfile(odir,sprintf('%s_%05d%s',n,i1,e));
            %     ni            = nifti;
            %     ni.dat        = file_array(Voo{i},V(i).dim(1:3),V(i).dt,0,V(i).pinfo(1),V(i).pinfo(2));
            %     ni.mat        = V(i).mat;
            %     ni.mat0       = V(i).mat;
            %     ni.descrip    = [V(i).descrip sprintf(' - %d',i)];
            %     create(ni);
            %     ni.dat(:,:,:) = V(i).private.dat(:,:,:,i); %spoooky
            %     overloading of normal assignment to write to disk on this
            %     line
            %spm_progress_bar('Set',i);
        end
        %spm_progress_bar('Clear');
        
        if nargout
            %     Vo = spm_vol(char(Voo));
            Vo = char(Voo);
        end
    end

    function val=readJson(fname)
        %someone on matlab answers wrote this
        fid = fopen(fname);
        raw = fread(fid,inf);
        str = char(raw');
        fclose(fid);
        val = jsondecode(str);
    end
end
