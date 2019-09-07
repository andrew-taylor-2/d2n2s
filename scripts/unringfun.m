function dwi_obj = unringfun(dwi_obj)
% operates on "objects" (structures)

%check dependency
if isempty(which('unring'))
    unringdir=fullfile(fileparts(fileparts(which(mfilename))),'utilities','external','unring_lib','matlab');
    if exists(fullfile(unringdir,['ringRm.' mexext]),'file')
        addpath(unringdir)
    elseif exists(exists(fullfile(unringdir,'unring.m'),'file'))
        error('need to compile ringRm.cpp into mex for unring dependency')
    else
        error('missing unring dependency')
    end
end


parfor j=1:length(dwi_obj)
    dwi_obj(j).img=unring(dwi_obj(j).img);
    dwi_obj(j).img(isnan(dwi_obj(j).img))=0;
    dwi_obj(j).hdr.dt=[64 0];
end
end
