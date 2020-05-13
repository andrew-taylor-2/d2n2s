function getslspec(in,outname)
% operates on either "objects" (structures) or filenames, depending on "in" variable
%if it's a field of ur a d2n2s structure (example:"in=dki(1).json"), just get it directly from there.
%if it's a .json file, do the FSL command
%
%variable input type, output is a file that can be given directly to FSL


if isstruct(in) %if the user gives a JSON structure
    %check if the manufacturer is siemens and then get the value straight
    %from the field
    if isfield(in,'SliceTiming')
        slicetimes=in.SliceTiming;
    else
        error('Though you''ve given this function a structure, there''s no field for SliceTiming.')
    end
    
    try manu=in.Manufacturer; catch; end
    
elseif isstring(in) || ischar(in) %if the user gives a filename
    in=char(in);
    
    %this next block is lifted straight from fsl's website
    fp = fopen(in,'r');
    fcont = fread(fp);
    fclose(fp);
    cfcont = char(fcont');
    i1 = strfind(cfcont,'SliceTiming');
    i2 = strfind(cfcont(i1:end),'[');
    i3 = strfind(cfcont((i1+i2):end),']');
    cslicetimes = cfcont((i1+i2+1):(i1+i2+i3-2));
    slicetimes = textscan(cslicetimes,'%f','Delimiter',',');
else 
    error('idk bro')
end
    
if exist('manu','var') && ~contains(manu,'Siemens'); warning('This doesn''t appear to be a Siemens-created .json file, which is an assumption of this script');end


%this bit is from fsl's website as well
[sortedslicetimes,sindx] = sort(slicetimes);
mb = length(sortedslicetimes)/(sum(diff(sortedslicetimes)~=0)+1);
slspec = reshape(sindx,[mb length(sindx)/mb])'-1;
dlmwrite(outname,slspec,'delimiter',' ','precision','%3d');
