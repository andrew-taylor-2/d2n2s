function flags = make_flags(read_write_str,varargin)

%make sure you know whether you're reading or writing
p_init=inputParser;
addRequired(p_init,'read_write_str', @(str) sum(strcmpi({'read','write'},str))==1 )
parse(p_init,read_write_str)
read_write_str=p_init.Results.read_write_str;

p=inputParser;

if strcmpi(read_write_str,'read')
    addOptional(p,'pick','',@(x) ischar(x)||isstring(x) )
    addOptional(p,'b0',1,@(x) isnumeric(x) && sum(x==[0 1])==1 )
    addOptional(p,'no','',@(x) ischar(x)||isstring(x) )
    %the 'no' option is going to become noload and nofns

    
    
elseif strcmpi(read_write_str,'write')
    
    addOptional(p,'dt',1, @(x) isequal(x,1) || isequal(x,-1) || (isequal(x(1),0) && sum(x(2)==[256 2 4 512 8 768 16 64])==1))
    addOptional(p,'gr',1,@(x) isnumeric(x) && sum(x==[0 1])==1 )
    addOptional(p,'vol',4,@(x) isnumeric(x) && sum(x==[3 4])==1 )
    addOptional(p,'nfn',0,@(x) isnumeric(x) && sum(x==[0 1])==1 )
    
end

parse(p,varargin{:})
flds=fieldnames(p.Results);
for i=1:numel(flds)
    flags.(flds{i})=p.Results.(flds{i});
end