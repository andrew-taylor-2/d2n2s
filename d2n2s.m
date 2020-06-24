function dwi=d2n2s(dcm2niixd_folder,varargin)
% operates on files/folders
%

%flags.pick - if this is specified, this nii filename will be chosen to have it's
%header/image information loaded. just in case, for instance, one has
%many nii files within their dcm2niixd_folder

%flags.glob - like flags.pick, but the string will be entered into dir and
%then fnify2 to pick the file. NOTE that while that while "pick"ing a .nii,
%users will enter the full path fn, users will "glob" only the name and
%extension of a file within the dcm2niixd_folder. The treatment of these
%two fields will hopefully be made more consistent in future updates.

%flags.no  - tells the program not to read certain things. could be used to
%speedup if you don't need that info. could be used to avoid errors an
%incorrect data loading. takes a single string that should contain 'bvec',
%'bval', 'img','nii', or 'json'. 'img' avoids loading the voxel data, 'nii'
%avoids loading the header and voxel data. For example: flags.no='bvecbvalnii' if
%all you want is json info

%flags.b0 - should we assume that folders not containing bvecs/bvals are
%b0s? default NO (I no longer use mostly diffusion images lel), but this could be changed by the contents of a json
%file, if present

%I found in making this function that there was a tradeoff between clarity
%about what the script was doing and making the script try its very best to
%give the user what they want. I chose to try to give the user what they
%want.

%it also might be useful to make the dir calls a wrapper, so that d2n2s can
%be made more general

% I need to make the new try <assign .fns field> statements their own
% function, it's making this function look bloated

%% random input sanitizing
% I SHOULD STICK FLAGS INTO make_flags() TO SANITIZE IT

%new line that should make things way more convenient...
flags=make_flags('read',varargin{:});

% if ~exist('flags','var')
%     flags=[];
% end


if isstring(dcm2niixd_folder)
    dcm2niixd_folder=char(dcm2niixd_folder);
end
if strcmp(dcm2niixd_folder(end),filesep) %strcmp for right now works on struct, let's see if that stays in later versions...
    dcm2niixd_folder=dcm2niixd_folder(1:end-1); %this isn't actually necessary since dir will ignore double separators, but I'm doing this for clarity
end

if isfield(flags,'pick') && ~isempty(flags.pick) && isstring(flags.pick)
    flags.pick=char(flags.pick);
end

if isfield(flags,'glob') && ~isempty(flags.glob) && isstring(flags.glob)
    flags.glob=char(flags.glob);
end

%don't have both pick and glob
if (isfield(flags,'pick') && ~isempty(flags.pick)) && (isfield(flags,'glob') && ~isempty(flags.glob))
    error('can''t have both pick and glob')
end



%these next two have to be initialized
if ~isfield(flags,'no') || isempty(flags.no) 
    flags.no='';
end

if ~isfield(flags,'b0') || isempty(flags.b0) 
    flags.b0=0;
end

if ~isfield(flags,'gz') || isempty(flags.gz) 
    flags.gz=0;
end


%check if they've picked or globbed a file,
using_pick=(isfield(flags,'pick') && ~isempty(flags.pick));
using_glob=(isfield(flags,'glob') && ~isempty(flags.glob));

%DETERMINE input type
%determine if first-position input is folder, file, or dir object
% if it's a dir object, make it into file or folder. if it's a file, put it
% into pick. if it's a folder, move it along as dcm2niixd_folder, to be
% used later as pp

% guessing dir will be easiest to tell or rule out
intype='';
if ~using_pick && ~using_glob
    
    % guessing DIR will be easiest to tell or rule out -- do that first
    if isstruct(dcm2niixd_folder) && isfield(dcm2niixd_folder,'name') && isfield(dcm2niixd_folder,'folder')
%         intype='dirobj';
        
        %we don't want to leave our input as a dir object, we want to
        %direct it into the file or folder stream. 
        
        temp=fnify2(dcm2niixd_folder);
        if exist(temp,'dir')
            intype='dir';
            dcm2niixd_folder=temp;
            
        elseif exist(temp,'file')
            intype='file';
            flags.pick=temp;
            using_pick=true;
            
        else
            error('your dir input isn''t a file or folder! might not exist')
        end
                
        
    %next try FILEname -- dir handling will be doing nothing as the rest of
    %the program already does dir handling
    elseif exist(dcm2niixd_folder,'file') && ~exist(dcm2niixd_folder,'dir')
        intype='file';
        flags.pick=dcm2niixd_folder;
        using_pick=true;
    else
        intype='dir';
    end
end




%get files with same name as pick or globbed file

if using_pick || using_glob
    
    if using_pick
        %look for files with the same "name"
        [pp,nn,ee]=fileparts(flags.pick);
    end
    
    if using_glob
        dir_out=dir(fullfile(dcm2niixd_folder,flags.glob));
        dir_out_fn=fnify2(dir_out);
        [pp,nn,ee]=fileparts(dir_out_fn);
    end
    
    %nn might contain
    if (length(nn)>4) && strcmp(nn(end-3:end),'.nii')
        nn=nn(1:end-4);
    end
    
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
    
    %yes, the logic flow is weird and gross looking. if only there were an
    %"if,else,finally"?
    
    %check gz, if so set gz to 1
    if contains(ee,'gz','IgnoreCase',true)
        flags.gz=1;
    end
    
end



%the next section was one of the first sections, the values got replaced
%later. I'm moving it later, but not assigning values if they already exist

%in earlier versions, "pp" was dcm2niixd_folder. This was easily replacable
%as long as i handled the case where we didn't pick or glob i.e. folder
%input 
if strcmp(intype,'dir')
    pp=dcm2niixd_folder;
end

if ~exist('dbvec','var')
    
    dbvec=dir([pp filesep '*.bvec']);
    dbvec2=dir([pp filesep '*.bvecs']);
    dbvec=[dbvec;dbvec2]; %just in case
end

if ~exist('dbval','var')
    dbval=dir([pp filesep '*.bval']);
    dbval2=dir([pp filesep '*.bvals']);
    dbval=[dbval;dbval2]; %just in case
end

if ~exist('djson','var')
    djson=dir([pp filesep '*.json']);
end









%% WARNINGS
if length(dbvec)>1
    if ~contains(flags.no,'bvec','IgnoreCase',1)
        warning('found more than one bvec file. using the ''first''.')
    end
    dbvec=dbvec(1);
end

if length(dbval)>1
    if ~contains(flags.no,'bval','IgnoreCase',1)
        warning('found more than one bval file. using the ''first''.')
    end
    dbval=dbval(1);
end
% if length(dnii)>1
%     warning('found more than one nii file. using the ''first''.')
% end
if length(djson)>1
    if ~contains(flags.no,'json','IgnoreCase',1)
        warning('found more than one json file. using the ''first''.')
    end
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
            warning(['your bvecs are square, so the orientation convention is ambiguous to this function. ' newline 'assuming they are in fsl/dcm2niix convention. user may need to transpose.'])
        else
            error('none of your bvec matrix dimensions are 3, unable to parse')
        end
        
        % ASSIGNMENT TO STRUCT
        [dwi(1:length(bvecs2)).bvec]=bvecs2{:};
    end
end

%% NIFTI FILE

%what if pick or glob fails? exit.
if (using_pick && ~exist(flags.pick,'file')) || (using_glob && isempty(dir_out))
    error('pick or glob failed, exiting')
end

% get nifti name with dir or flags.pick
if ~using_pick && ~using_glob %if the user has not opted to specify the nii fn himself
    
    dnii=dir([pp filesep '*.nii']);
    
    if isempty(dnii)
        dnii=dir([pp filesep '*.nii.gz']);
        if ~isempty(dnii)
            warning('couldn''t find a .nii file, but found a .nii.gz -- using')
            flags.gz=1;
        end
    end
    
    
elseif using_pick % if the user HAS specified the nii fn himself
    dnii=dir(flags.pick);
    
elseif using_glob
    dnii=dir_out;
end

nii_file=fnify2(dnii);
if ~isempty(nii_file) && strcmp('.gz',nii_file(end-2:end)) && flags.gz==1
    do=get_anonymous_functions;
    niigz_file=nii_file;
    nii_file=do.gunzip_and_rename_no_delete(nii_file);
    nii_file=do.move_and_rename(nii_file,[tempdir choose_output(@() fileparts(nii_file),2) dicomuid '.nii']);
end
    
    
    
if ~isempty(nii_file)
    
    % if the user doesn't want the image data NOR the header, don't do
    % anything in this section
    if ~contains(flags.no,'nii','IgnoreCase',1)
        
        %load hdr
        hdr=spm_vol(nii_file);
        
        % sanitize
        if exist('dwi','var')
            if length(hdr)~=length(dwi)
                warning(['spm_vol thinks your 4D volume has a different number of images from your bvals and/or bvecs file.' newline 'this function might not behave as intended'])
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
            
            if numel(dwi(1).hdr.pinfo)>=2 && any(dwi(1).hdr.pinfo(1:2)~=[1;0])
                warning('the .hdr.pinfo field indicates that the "true" voxel values might be different from the raw values stored (e.g. the values might be scaled). One consequence of this is that values written to the image (after using d2n2s_write or spm_write_vol) may be only approximations of those that were assigned.')
            end
            
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
            % it was actually important somehow. I'm still not sure exactly why this is, but
            % I've checked that this is functioning well empirically and it
            % is so.... it stays.
            
            % new new comment: 
            % this line does seem to alert spm to the fact that it doesn't
            % support .nii.gz files and causes an error but... they're "not
            % supported" anyway so whatever
            
            % new new new comment:
            % I'm commenting it. It doesn't make sense to use a hack I
            % don't understand to reduce the functionality of this program
            
            % shiny new comment:
            % As earlier comments have noted, commenting allows you to use
            % .nii.gz files. However, .nii.gz files will inexplicably get a
            % obj.hdr.dat field (((possibly because the contrived .fn that is
            % handed to .private.dat.fname is always .nii and not .nii.gz,
            % and thus escapes the "unknown filename extension" error in
            % line 404 of nifti/subsasgn>assigndat (wait, you would think the hack would help files slip by... what) ))), while normal .nii
            % won't have this (both types do get a obj.hdr.private.dat
            % field, however). This will cause a "subscripted assignment
            % between dissimilar structures" error in spm_reslice_at.m,
            % specifically at 
            
            % for i=1:numel(P2)
            % P(i)=P2(i).hdr; %field 'hdr' must exist
            % end

            % it seems like this could easily be remedied, but I have yet
            % again been spooked out and reminded that "compressed nifti
            % files are not supported" and the hassle of converting .nii.gz
            % files is less than the stress of always wondering if my code
            % is breaking silently but deadlyily.
            
            dwi(i).hdr.private.dat.fname=v0(i,:); %old old comment: this line is included to protect the user, as assigning to dat(:) at all in matlab will overwrite whatever is contained in hdr.private.dat.fname. this is a property of the @file_array object. and it's spooky.
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

%assign niigz if you have gone through the gz pipeline -- it's nice to have
%a path to the original data
if exist('niigz_file','var')
try
    for i=1:length(dwi)
        dwi(i).fns.niigz=niigz_file;
    end
catch
    warning('failed assigning niigz fns, continuing')
end
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
    warnings=[false false];
    if ~contains(flags.no,'bvec','IgnoreCase',1)
        
        % assign to struct
        [dwi.bvec]=deal([0;0;0]);
        
        %warn
        warnings(1)=true;
    end
    
    % if we want bvals filled in
    if ~contains(flags.no,'bval','IgnoreCase',1)
        
        % assign to struct
        [dwi.bval]=deal(0);
        
        %warn
        warnings(2)=true;
    end
    
    % warnings -- code could be made prettier or cooler but I think it
    % would mean slowing down and/or making the warning message harder to
    % read
    if contains(dcm2niixd_folder,'b0','IgnoreCase',1)
        warning_prefix=['there are no bvals or bvecs, and the folder you''ve given this program contains the string ''b0''' newline];
    else
        warning_prefix=['there are no bvals or bvecs' newline 'ALSO the folder you have supplied doesn''t have ''b0'' in the name. Still,' newline];
    end
    
    warning_suffix=[' if this isn''t desired, set flags.b0 to [0]'];
    if all(warnings)
        warning([warning_prefix 'setting bvals and bvecs to 0.' warning_suffix])
    elseif warnings(1)
        warning([warning_prefix 'setting bvecs to 0.' warning_suffix])
    elseif warnings(2)
        warning([warning_prefix 'setting bvals to 0.' warning_suffix])
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

function do = get_anonymous_functions

%inline conditional
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();

%index even unnamed expressions
do.paren=@(x,varargin) x(varargin{:});
do.curly=@(x,varargin) x{varargin{:}};

%add onto a filename -- these have been tested, and they're about 3 times
%slower than their multiline function counterparts. This is a negligible
%portion of the total runtime of this routine, though. I'm using these
%anonymous functions in the first place for clarity and readability of
%code. I think they're worth it
do.append=@(x,str) fullfile(choose_output(@() fileparts(x),1),[choose_output(@() fileparts(x),2) str choose_output(@() fileparts(x),3)]);
do.prepend=@(x,str) fullfile(choose_output(@() fileparts(x),1),[str choose_output(@() fileparts(x),2) choose_output(@() fileparts(x),3)]);

%do something and return the out name
do.move_and_rename=@(in,out) do.curly({movefile(in,out),out},2);
do.copy_and_rename=@(in,out) do.curly({copyfile(in,out),out},2);
do.gunzip_and_rename=@(in) iif( ~ischar(in) || numel(in)<4,  @() do.curly({in, warning('invalid gunzip input,returning in unchanged')},1), ... %warning for invalid inputs
                                strcmp(in(end-2:end),'.gz'), @() do.curly({gunzip(in),void(@() delete(in)),in(1:end-3)},3), ... % if it's named ".gz", gunzip and return out name
                                true,                        @() in); % if it's not just return the in name

do.gunzip_and_rename_no_delete=@(in) iif( ~ischar(in) || numel(in)<4,  @() do.curly({in, warning('invalid gunzip input,returning in unchanged')},1), ... %warning for invalid inputs
                                strcmp(in(end-2:end),'.gz'), @() do.curly({gunzip(in),in(1:end-3)},2), ... % if it's named ".gz", gunzip and return out name
                                true,                        @() in); % if it's not just return the in name

function o=void(f)
o=1;
f();
end
end
end
