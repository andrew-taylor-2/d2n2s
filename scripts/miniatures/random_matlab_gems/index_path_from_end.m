function out=index_path_from_end(path1,index)
curly=@(x,varargin) x{varargin{:}};
cell_path=curly(textscan(path1,'%s','Delimiter','/'),1);
out=fullfile(cell_path{end+1-flip(index)});
end
