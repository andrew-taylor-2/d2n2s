function flags = make_flags(read_write_str,varargin)

%make sure you know whether you're reading or writing
p_init=inputParser;
addRequired(p_init,'read_write_str', @(str) sum(strcmpi({'read','write'},str))==1 )
parse(p_init,read_write_str)
read_write_str=p_init.Results.read_write_str;
clear p_init
p=inputParser;

%annoyingly, I have to do this again or it'll mess up
addRequired(p,'read_write_str', @(str) sum(strcmpi({'read','write'},str))==1 )
%in this weird case it's actually useful for them not to programmatically
%assign things to my variables....

if strcmpi(read_write_str,'read')
    addParameter(p,'pick','',@(x) ischar(x)||isstring(x) )
    addParameter(p,'b0',1,@(x) isnumeric(x) && sum(x==[0 1])==1 )
    addParameter(p,'no','',@(x) ischar(x)||isstring(x) )
    %the 'no' option is going to become noload and nofns

    
    
elseif strcmpi(read_write_str,'write')
    
    addParameter(p,'dt',1, @(x) isequal(x,1) || isequal(x,-1) || (isequal(x(1),0) && sum(x(2)==[256 2 4 512 8 768 16 64])==1))
    addParameter(p,'gr',1,@(x) isnumeric(x) && sum(x==[0 1])==1 )
    addParameter(p,'vol',4,@(x) isnumeric(x) && sum(x==[3 4])==1 )
    addParameter(p,'nfn',0,@(x) isnumeric(x) && sum(x==[0 1])==1 )
    
end

parse(p,read_write_str,varargin{:})
flds=fieldnames(p.Results);
for i=1:numel(flds)
    if ~strcmpi(flds{i},'read_write_str') %this is a meta variable, don't want it in the output
        flags.(flds{i})=p.Results.(flds{i});
    end
end