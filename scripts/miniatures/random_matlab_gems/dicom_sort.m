function dicom_sort(fn_o)
%base_folder = uigetdir('', 'Input folder:');

%some of these lines (and the inspiration) were written by Russell Glenn

subject_folder = fn_o;

% move all files
%     eval('!for /f "delims=" %f in (''dir/s/w/b/A:D'') do move "%f\*.*"')

%     % above is the old command that called the windows command line. it
%     grabs all files inside the specified folder (even nested in folders)
%     and moves them to the subject folder

files=clean_dir2(dir(fullfile(subject_folder,'/**/*')));
files([files.isdir])=[]; %don't move any directories -- they'll be deleted later

for j=1:numel(files)
    try
        movefile(fnify2(files(j)),subject_folder);
    catch %this will catch when the dicom that is grabbed is already in the current folder
    end
end

% delete folders
%     eval('!for /f "delims=" %f in (''dir/w/b/A:D'') do rd/q/s "%f"')
%     
%     % above is the old windows command line. it deleted all nested
%     folders. This code contains an untested section that changes
%     "clean_dir" to remove only . and .. from the dir objects so you
%     wouldn't have to kludge to try to delete ds store and any other file
%     names starting with a dot wouldn't cause a problem.

emptydirs=clean_dir2(dir(fullfile(subject_folder,'/**/')));
emptydirs(~[emptydirs.isdir])=[];

%     for j=numel(emptydirs):-1:1; try rmdir(fnify(emptydirs(j))); catch;end;end
%     nuisance_dir=dir(fullfile(subject_folder,'/**/.DS_Store'));
%     for j=numel(nuisance_dir); delete(fnify(nuisance_dir(j)));end

%remove empty dirs
for j=numel(emptydirs):-1:1 %don't think this actually needs to go backwards
    rmdir(fnify2(emptydirs(j)));
end


% get directory listing
filelist = clean_dir(dir(subject_folder));
filelist = {filelist.name}';

for j=length(filelist):-1:1 %goes backwards so that deleting won't mess up indexing
    if contains(filelist{j},'.bvec') || contains(filelist{j},'.bval') || contains(filelist{j},'.json') || contains(filelist{j},'.nii') || contains(filelist{j},'.info') || contains(filelist{j},'.nfo')% this leaves any previously converted files in the subject/given folder and doesn't sort them. They also shouldn't cause an error. A better way to do this (albeit probably slower) would be to check the file extension against a list of typical file extensions. If it doesn't have a file extension or does have a .dcm ending, don't remove it from filelist
        filelist(j)=[];
    end
end

serieslist = {};

% move file to folder named after 'series'
ids=0;
for j = 1:length(filelist)
    
    fn = fullfile(subject_folder, filelist{j});
    
    info = dicominfo(fn);
    
    if isfield(info, 'SeriesDescription')
        
        SeriesID = [info.SeriesDescription '_' num2str(info.SeriesNumber)];
        
    else
        try
            SeriesID = num2str(info.SeriesNumber);
        catch
            ids=ids+1;
            SeriesID = ['noID' num2str(ids)]
        end
        
    end
    
    target_folder = fullfile(subject_folder, SeriesID);
    
    %only make a new folder for a series if you need to
    if isempty(strcmpi(SeriesID, serieslist)) || ~any(strcmp(SeriesID, serieslist))
        
%         fprintf('Found new series %s for subject %s\n', SeriesID, subject_list{i})
        
        mkdir(target_folder);
        
        serieslist = [serieslist(:)' SeriesID];
        
    end
    
    movefile(fn, target_folder);
    
end


end
function fn = fnify2(dir_object,bhvr_bin)
% bhvr_bin = 1: error when dir_object is empty
% bhvr_bin = 0: fn is empty string when dir_object is empty
if ~exist('bhvr_bin','var')
    bhvr_bin=0;
end
if isempty(dir_object) && bhvr_bin
    error('your dir call was empty, so I can''t make a filename from it')
elseif isempty(dir_object) && ~bhvr_bin
    fn='';
    return
elseif numel(dir_object)>1
    warning(['numel of object from dir call is greater than 1.\n -- it''s ' num2str(numel(dir_object)) ' -- using the least recently modified file.']);
end
dir_object=dir_object(choose_output(@() min(cat(1,dir_object.datenum)),2)); 

fn=[dir_object.folder filesep dir_object.name];
end

function cdir=clean_dir(dir) %gets rid of current folder, parent folder, and hidden files.
bdir=arrayfun(@(x) x.name(1)=='.',dir);
cdir=dir(~bdir);
end

function cdir=clean_dir2(dir) %only gets rid of current folder and parent folder
bdir=arrayfun(@(x) strcmp(x.name,'.') || strcmp(x.name,'..') ,dir);
cdir=dir(~bdir);
end
