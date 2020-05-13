function out=index_path2(path1,index)
curly=@(x,varargin) x{varargin{:}};

cell_path=curly(textscan(path1,'%s','Delimiter',filesep),1);
out=fullfile(cell_path{index+1});

%check for annoying case:
if ~exist(out,'file') && exist([filesep out],'file')
    warning('watch out indexing path, you might be leaving out the root folder')
end

end

%two issues: 

%1. the way this is written, '/', the root folder, isn't
%considered to be a folder. It's not counted in the indexing, and, if
%indexing from the beginning, it's not included. That last bit is a tricky
%problem. just cost me like 5-10 minutes scratching my head. 

%2. we don't know enough about the use case. The user could be indexing an
%absolute or full path, and the user could be putting together a file that
%does exist or doesn't yet exist. Instead of writing this to be right, I'm
%just going to include a warning for a common pitfall and keep the function
%how I'm familiar with.


% below is another version of the script. I scrapped it because it only
% works well for absolute paths

%this function used to behave like the ~isunix case always. but it can get
%confusing when you don't consider the root folder to be a folder to
%index.... I ended up writing tons of files to `pwd`/Users/<myuser>/....
%because I didn't include a root folder...


% curly=@(x,varargin) x{varargin{:}};
% if isunix
%     cell_path=curly(textscan(path1,'%s','Delimiter',filesep),1);
%     cell_path{1}='/';
%     out=fullfile(cell_path{index});
% else
%     cell_path=curly(textscan(path1,'%s','Delimiter',filesep),1);
%     out=fullfile(cell_path{index+1});
% end
% 
% %this function used to behave like the ~isunix case always. but it can get
% %confusing when you don't consider the root folder to be a folder to
% %index.... I ended up writing tons of files to `pwd`/Users/<myuser>/....
% %because I didn't include a root folder...




