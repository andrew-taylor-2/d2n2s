function dwi=d2n2s(dcm2niixd_folder,flags)
% operates on files/folders
%
%flags.pick - if this is specified, this nii filename will be chosen to have it's
%header/image information loaded. just in case, for instance, one has
%many nii files within their dcm2niixd_folder

%flags.no  - tells the program not to read certain things. could be used to
%speedup if you don't need that info. could be used to avoid errors an
%incorrect data loading. takes a single string that should contain 'bvec',
%'bval', 'img','nii', or 'json'

%I found in making this function that there was a tradeoff between clarity
%about what the script was doing and making the script try its very best to
%give the user what they want. I chose to try to give the user what they
%want.

%it also might be useful to make the dir calls a wrapper, so that d2n2s can
%be made more general

% I need to make the new try <assign .fns field> statements their own
% function, it's making this function look bloated

%% random input sanitizing

if ~exist('flags','var')
    flags=[];
end

if isstring(dcm2niixd_folder)
    dcm2niixd_folder=char(dcm2niixd_folder);
end
if isfield(flags,'pick') && ~isempty(flags.pick) && isstring(flags.pick)
    flags.pick=char(flags.pick);
end

%these next two have to be initialized 
if ~isfield(flags,'no') || isempty(flags.no) 
    flags.no='';
end

if ~isfield(flags,'b0') || isempty(flags.b0) 
    flags.b0=1;
end


%% GRAB FILES

dbvec=dir([dcm2niixd_folder filesep '*.bvec']);
dbvec2=dir([dcm2niixd_folder filesep '*.bvecs']);
dbvec=[dbvec;dbvec2]; %just in case
dbval=dir([dcm2niixd_folder filesep '*.bval']);
dbval2=dir([dcm2niixd_folder filesep '*.bvals']);
dbval=[dbval;dbval2]; %just in case
%dnii=dir([dcm2niixd_folder filesep '*.nii']); %dnii gets set later and is part of a conditional
djson=dir([dcm2niixd_folder filesep '*.json']);

%if they've picked a file,
if isfield(flags,'pick') && ~isempty(flags.pick) 
    %look for files with the same "name"
    [pp,nn,ee]=fileparts(flags.pick);
    
    matching_bvec_name=[pp filesep nn '.bvec'];
    if exist(matching_bvec_name,'file')
        dbvec=dir(matching_bvec_name);
    end
       
    matching_bval_name=[pp filesep nn '.bval'];
    if exist(matching_bval_name,'file')
        dbval=dir(matching_bval_name);
    end
    
    matching_json_name=[pp filesep nn '.json'];
    if exist(matching_json_name,'file')
        djson=dir(matching_json_name);
    end
    
end
    
    

%% WARNINGS
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


%% BVALS
% prepare bvals and grab bval_file

if ~isempty(dbval)
    
    %grab file -- to be assigned to all dwis later
    bval_file=fnify2(dbval);
    
    if ~contains(flags.no,'bval','IgnoreCase',1)
        
        %load bvals
        bvals=load(bval_file);
        bvals2=num2cell(bvals);
        
        % ASSIGNMENT TO STRUCT
        [dwi(1:length(bvals2)).bval]=bvals2{:};
    end
end

%% BVECS 
% prepare bvecs and grab bvec_file

if ~isempty(dbvec)
    
    %grab file -- to be assigned to all dwis later
    bvec_file=fnify2(dbvec); % I can grab this here regardless of contains(flags.no,'bvec') bc there is just about 0 chance this would cause an error. potential assignment to fn can go at the end
    
    if ~contains(flags.no,'bvec','IgnoreCase',1)
        
        %load and sanitize bvecs
        bvecs=load(bvec_file);
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
        
        % ASSIGNMENT TO STRUCT
        [dwi(1:length(bvecs2)).bvec]=bvecs2{:};
    end
end

%% NIFTI FILE

% get nifti name with dir or flags.pick
if ~isfield(flags,'pick') || isempty(flags.pick) || ~exist(flags.pick,'file') %if the user has not opted to specify the nii fn himself
    dnii=dir([dcm2niixd_folder filesep '*.nii']);
    
    if isempty(dnii)
        dnii=dir([dcm2niixd_folder filesep '*.nii.gz']);
        warning('couldn''t find a .nii file, but found a .nii.gz -- using')
    end
    
    
elseif isfield(flags,'pick') && ~isempty(flags.pick) && ~isempty(dir(flags.pick)) % if the user HAS specified the nii fn himself
    dnii=dir(flags.pick);
end


if ~isempty(dnii)
    
    %grab file. will be put in dwi.fns.nii later
    nii_file=fnify2(dnii);
    
    % if the user doesn't want the image data NOR the header, don't do
    % anything in this section
    if ~contains(flags.no,'nii','IgnoreCase',1)
        
        %load hdr
        hdr=spm_vol(nii_file);
        
        % sanitize
        if exist('dwi','var')
            if length(hdr)~=length(dwi)
                warning('spm_vol thinks your 4D volume has a different number of images from your bvals file. this function might not behave as intended')
            end
        end
        
        %prep to assign header
        hdr2=cell(1,length(hdr));
        for i=1:length(hdr)
            hdr2{i}=hdr(i);
        end
        
        % ASSIGNMENT TO STRUCT
        [dwi(1:length(hdr2)).hdr]=hdr2{:};
        
        % skip actually reading of image if user wishes
        if ~contains(flags.no,'img','IgnoreCase',1)
            
            % read in images, prep to assign
            img=spm_read_vols(hdr);
            img2=num2cell(img,[1 2 3]);
            
            % ASSIGNMENT TO STRUCT
            [dwi.img]=img2{:};
        end
        
        %this function is kind of obsolete
        v0=spm_file_split2(nii_file);
        for i=1:size(v0,1)
            dwi(i).fn=v0(i,:);

            
            % old comment:
            % I don't remember exactly how the line below works anymore.
            % I'm keeping it just because it didn't hurt anything and I've
            % really thoroughly tested this function at this point...
            % Honestly I think p...d...fname gets overwritten by spm
            % anyway. If anyone is curiously reading this, ask me and I'll
            % tell you about it to the best of my understanding
            
            % new comment:
            % it was actually important somehow. I commented it, and later
            % on, when I tried to concatenate spm_vol objects which
            % previously had different datatypes and write, it messed up
            % and wrote noise. I'm still not sure exactly why this is, but
            % I've checked that this is functioning well empirically and it
            % is so.... it stays.
            
            dwi(i).hdr.private.dat.fname=v0(i,:); %this line is included to protect the user, as assigning to dat(:) at all in matlab will overwrite whatever is contained in hdr.private.dat.fname. this is a property of the @file_array object. and it's spooky.
        end
    end
end


%% JSON

if ~isempty(djson)
    
    %grab fn
    jsfn=fnify2(djson);
    
    if ~contains(flags.no,'json','IgnoreCase',1)
        
        % load vals
        jsvals=readJson(jsfn);
        % assignment to struct
        [dwi.json]=deal(jsvals); % the json values are assigned to all dwis, like .fns.blah 
    end
end


%% assign all fns (goes here bc each element of array has the same value)

%NII -- we're always going to have a nii fn, so never exclude this
%assignment
try
    for i=1:length(dwi)
        dwi(i).fns.nii=nii_file;
    end
catch
    warning('failed assigning nii fns, continuing')
end

if ~contains(flags.no,'fn','IgnoreCase',1)
    
    %BVEC
    if exist('bvec_file','var')
    try
        for i=1:length(dwi)
            dwi(i).fns.bvec=bvec_file;
        end
    catch
        warning('failed assigning bvec fns, continuing')
    end
    end
    
    %BVAL
    if exist('bval_file','var')
    try
        for i=1:length(dwi)
            dwi(i).fns.bval=bval_file;
        end
    catch
        warning('failed assigning bval fns, continuing')
    end
    end
    
    %JSON
    if exist('jsfn','var')
    try
        for i=1:length(dwi)
            dwi(i).fns.json=jsfn;
        end
    catch
        warning('failed assigning json fns, continuing')
    end
    end
end


%% ASSIGN METADATA when d2n2s doesn't provide it to avoid failure or manual assignment down the road

%some b0 series do not have bvals or bvecs -- for very specific cases,
%include these
if ~isequal(flags.b0,0) && exist('dwi','var') && isempty(dbvec) && isempty(dbval) % it's only really suspicious if both are missing.
    
    % a t1 would get to this point, though, so, if we have the resources, let's see if we can tell
    % this is a t1 a little better 
    try
    if isfield(dwi,'json')
        
        conds(1)=contains(dwi(1).json.SeriesDescription,'t1','IgnoreCase',1);
        conds(2)=contains(dwi(1).json.SeriesDescription,'mprage','IgnoreCase',1);
        
        conds(3)=contains(dwi(1).json.ProtocolName,'t1','IgnoreCase',1);
        conds(4)=contains(dwi(1).json.ProtocolName,'mprage','IgnoreCase',1);
        
        conds(5)=~any(strcmpi(dwi(1).json.ImageType,'diffusion'));

        if any(conds); return;end % if any of these conditions are true, "return" -- this will prevent assignement of b0 metadata
    end
    catch
    end
    
    % if we want bvecs filled in
    if ~contains(flags.no,'bvec','IgnoreCase',1)
        % assign to struct
        [dwi.bvec]=deal([0;0;0]);
    end
    
    % if we want bvals filled in
    if ~contains(flags.no,'bval','IgnoreCase',1)
        % assign to struct
        [dwi.bval]=deal(0);
    end
    
    % warning
    if contains(dcm2niixd_folder,'b0','IgnoreCase',1)
        warning('there are no bvals or bvecs, and the folder you''ve given this program contains the string ''b0''. setting bvals and bvecs to 0.')
    else
        warning('there are no bvals or bvecs\nALSO the folder you have supplied doesn''t have ''b0'' in the name. \nStill, setting bvals and bvecs to 0 as a best guess.')
    end
end


%% END OF MAIN

%% FUNCTIONS USED

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
