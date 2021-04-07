function out_fn=d2n2s_write(dwi,folder,name,varargin)
%take matlab dwi object and write it to file
% input is "object" (structure created by d2n2s), output is files
%from andrew's OOP scripts
%flags
%  flags.dt - decides behavior when combining data previously stored as
%      multiple different datatypes. 1 (default) saves all output as same dt as
%      most precise input dt, -1 saves all output as same dt as least precise
%      input dt, [0 n] or [0;n] where n is an type identifying number from
%      spm_type() ('types') saves outputs as this specified type.
%
%  flags.gr - 1 for writing gradient text files the cbi way (tall). default is 1. 0 for the fsl way
%  flags.vol - 3 for 3d files. 4 for a 4d file
%  flags.nfn - 1 for writing 3d volumes according to 'folder','name'.
%            - 0 for using the 'dwi.fn' field to write. 1 is default
%  flags.del - 1 for deleting the file you're about to write to. default 0

%if you wanted to use this on only one field, you could make another object
%with just the field. or use a segment of the object

% note: this is an SPM issue, but if you try to write to a folder and use
% the home folder '~' notation, e.g. folder='~/re/test'; spm is going to
% fail. just replace '~' with your home directory

%% handle inputs
flags=make_flags('write',varargin{:});

%
% if ~exist('flags','var')
%     flags=[];
% end

if ~isempty(folder) && ~exist(folder,'dir')
    try
        mkdir(folder)
    catch
    end
end

if ~isfield(flags,'del') || isempty(flags.del)
    flags.del=0;
end

%strip nii if necessary
if contains(name,'.nii') && strcmp(name(end-3:end),'.nii')
    name=name(1:end-4);
end

%% bvals
% if isfield(dwi,'bval') && ~isempty([dwi.bval])
if ~isempty([dwi.bval])
    bvalw=[dwi(:).bval];
    if numel(bvalw)<numel(dwi);
        warning('you might have some empty bvals; proceeding anyway')
    end
    fn=fullfile(folder,[name '.bval']);
    fid=fopen(fn,'w');
    fprintf(fid,'%u ',bvalw);
    fclose(fid);
end

%% bvecs
% if isfield(dwi,'bvec') && ~isempty([dwi.bvec])
if ~isempty([dwi.bvec])
    bvecs=[dwi(:).bvec];
    if length(bvecs)<numel(dwi)
        warning('some of your bvecs might be empty; proceeding anyway')
    end
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


%% nii file
% if isfield(dwi,'hdr') && isfield(dwi,'img')
fn=fullfile(folder,[name '.nii']);
if ~isfield(flags,'vol') || isempty(flags.vol) || ~any(flags.vol==[3,4])
    flags.vol=4;
end

%% if you're writing a 4d volume
if flags.vol==4
    fn=fullfile(folder,[name '.nii']); %this is named according to the function's inputs
    if exist(fn,'file') && flags.del~=1; warning('you''re writing to an image file that already exists; note that if you overwrite an existing image with an image containing less volumes, the extra volumes at the end will stay');end
    
    %% handle variable input dts
    
    %write to single out dt
    %note that none of this handles dt(2) ie big/little
    %endian-ness... I've never encountered a problem with this. if
    %you were to, you could just set
    %obj(:).hdr.dt(2)=spm_platform('bigend')
    
    %set flag if it hasn't been
    if ~isfield(flags,'dt') || isempty(flags.dt) || ~any(flags.dt(1)==[0 1 -1]) %if flag isn't set of is invalid
        %then make it 1
        flags.dt=1;
    end
    
    % spm_type types ordered by increasingly more precise, then increasing unsignedness (there has to be a better way to say that).
    types=[256 2 4 512 8 768 16 64];
    
    if isequal(flags.dt(1),0) %if the user has specified this option, you can skip all of the input datatype checking and handling
        %use whatever datatype the user specifies
        assert(any(flags.dt(2)==types),'specified type has to be valid spm_type')
        type_to_write=flags.dt(2);
        
    else
        
        %if the user hasn't handpicked an output datatype, then you
        %have to handle the input datatypes and decide which to use
        %between them.
        
        %get all the dts
        try dts=arrayfun(@(x) x.hdr.dt(1),dwi);catch; error('all .hdr.dt(1) elements must exist and be the same type');end
        assert(isnumeric(dts),'all .hdr.dt(1) elements must be numeric')
        
        %get indices
        [type_index_of_dt,~]=ind2sub(size((dts(:)==types)'),find((dts(:)==types)')); %alternative is to use find on vectors using arrayfun or cells or something
        assert(numel(type_index_of_dt)==numel(dts),'all .hdr.dt(1) elements must be in the set of spm types found in spm_type')
        
        % homogenize datatypes based on value in flags.dt
        if isequal(flags.dt(1),1) % if flags.dt==1 or it hasn't been set properly
            %then set all dt to the most precise in inputs
            type_to_write=types(max(type_index_of_dt));
            
        elseif isequal(flags.dt(1),-1) %if flags.dt==-1
            %then set all dt to the least precise in inputs
            type_to_write=types(min(type_index_of_dt));
        end
    end
    
    %set all obj.hdr.dt(1) to this datatype
    
    %% set values and write 4d
    if flags.del==1
        warning('off','MATLAB:DELETE:FileNotFound')
        delete(fn)
        warning('on','MATLAB:DELETE:FileNotFound')
    end
    for i1=1:length(dwi)
        dwi(i1).hdr.fname=fn;
        dwi(i1).hdr.n=[i1,1];
        dwi(i1).hdr.dt(1)=type_to_write;
        
        spm_write_vol(dwi(i1).hdr,dwi(i1).img);
    end
end
%% if you're writing 3d files
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
            
            %get the most succinct name
            if length(dwi)>=100
                formstring='03';
            elseif length(dwi)>=10
                formstring='02';
            else
                formstring='01';
            end
            
            
            fn=fullfile(folder,[name sprintf(['_%' formstring 'd.nii'],i2)]); %the volumes are programmatically named, similar to how
            %spm_split_file would name them
            dwi(i2).hdr.fname=fn;
            dwi(i2).hdr.n=[1,1];
            
            spm_write_vol(dwi(i2).hdr,dwi(i2).img);
        end
    end
end
% end

%make the output
if flags.vol==3
    %     namescell=arrayfun(@(x) x.hdr.fname,dwi,'un',0);
    %     for i=1:numel(namescell{1}) %assume they're all the same length bc that isss how it works
    %         temp=cellfun(@(x) x(i),namescell,'un',0);
    %         samemask(1,i)=isequal(temp{:});
    %     end
    %     out_fn=
    
    %that was the beginning of my solution to the "make a pattern from the
    %filenames and assign as output" problem, but then I remembered I've
    %already done this before
    out_fn=make_a_pattern(arrayfun(@(x) x.hdr.fname,dwi,'un',0));
else
    out_fn=fn;
end
%% json file
% if isfield(dwi,'json')
empty_jsons=arrayfun(@(x) isempty(x.json),dwi);
nonempty_jsons_inds=find(~empty_jsons);
if ~all(empty_jsons)
    try; if numel(dwi)>1 && ~isequal(dwi.json)
            warning('Not every json in the dataset to be written is equal to one another. Writing only the first nonzero json element''s .json field')
        end; catch; end
    jsontext=jsonencode(dwi(nonempty_jsons_inds(1)).json);
    fn=fullfile(folder,[name '.json']);
    fid=fopen(fn,'w');
    fprintf(fid,'%s',jsontext);
    fclose(fid);
end
% end
end

function the_pattern=make_a_pattern(patterned_cellstr)
%makes a pattern by matching characters forward into the string until a
%divergence, then backward, then joins the two with an asterisk.
%
%takes in a cell array of strings

%takes about 0.002 seconds for numel(patterned_cellstr)==10 and time increases
%linearly with its length and with the length of the beginning and ending
%strings

max_identicality_index=-1; %start 1 away from lowest possible ("doesn't end up indexing anything") bc it gets incremented once irrespective of the data
samee=true;
while samee==true
    max_identicality_index=max_identicality_index+1;
    letters=cellfun(@(x) x(max_identicality_index+1),patterned_cellstr,'un',0);
    samee=all(isequal(letters{:}));
end

min_identicality_index=2; %start 1 away from highest possible ("doesn't end up indexing anything") bc it gets incremented once irrespective of the data
samee=true;
while samee==true
    min_identicality_index=min_identicality_index-1;
    letters=cellfun(@(x) x(numel(x)+min_identicality_index-1),patterned_cellstr,'un',0);
    samee=all(isequal(letters{:}));
end

the_pattern=[patterned_cellstr{1}(1:max_identicality_index) '*' patterned_cellstr{1}(end+min_identicality_index:end)];
end
