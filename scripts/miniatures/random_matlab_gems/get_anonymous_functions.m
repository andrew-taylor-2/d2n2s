
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
do.gzip_and_rename=@(in) iif( ~ischar(in) || numel(in)<4,  @() do.curly({in, warning('invalid gunzip input,returning in unchanged')},1), ... %warning for invalid inputs
                                ~strcmp(in(end-2:end),'.gz'), @() do.curly({gzip(in),void(@() delete(in)),[in '.gz']},3), ... % if it's named ".gz", gunzip and return out name
                                true,                        @() in); % if it's already a gz just return the in name
function o=void(f)
o=1;
f();
end
end