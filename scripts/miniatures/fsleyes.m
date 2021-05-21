function fsleyes(DWI_obj)
%simple, effective. Probably.
% system(['fsleyes ' DWI_obj(1).fns.nii])
if length(DWI_obj)==1
    cmd=['fsleyes ' DWI_obj(1).fns.nii];
else
    fn_cellstr=arrayfun(@(x) x.fns.nii,DWI_obj);
    cmd=['fsleyes ' sprintf('%s ',string(fn_cellstr))];
    disp(cmd)
end
system(cmd)
end