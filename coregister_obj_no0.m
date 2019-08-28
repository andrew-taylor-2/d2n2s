
function [new_source,M]=coregister_obj_no0(target_object_seg,source_obj_seg)
%target object seg should be a dwi object with one element. 
if length(target_object_seg)>1; error('unintended usage');end

% coregistration and reslicing parameters
estflg.cost_fun = 'nmi';
estflg.sep      = [4 2];
estflg.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
estflg.fwhm     = [7 7];
wrtflg        = spm_get_defaults('realign.write');
wrtflg.interp   = 1;
wrtflg.which    = [1 0];

target_object_seg.hdr.img=target_object_seg.img;
for i=1:length(source_obj_seg)
    source_obj_seg(i).hdr.img=source_obj_seg(i).img;
end

x  = spm_coreg_at_no0(target_object_seg.hdr, source_obj_seg(1).hdr, estflg); 

%unset .hdr.img for clarity 
target_object_seg.hdr=rmfield(target_object_seg.hdr,'img');
for i=1:length(source_obj_seg)
    source_obj_seg(i).hdr=rmfield(source_obj_seg(i).hdr,'img');
end


M  = inv(spm_matrix(x));
MM = zeros(4, 4, length(source_obj_seg));

%need to do my own spm_get_space
for j=1:length(source_obj_seg)
    %MM(:,:,j) = spm_get_space(deblank(fn_other(j,:)));
    MM(:,:,j)=source_obj_seg(j).hdr.mat;
end
for j = 1:length(source_obj_seg)
    source_obj_seg(j).hdr.mat=M*MM(:,:,j);
    source_obj_seg(j).hdr.private.mat=source_obj_seg(j).hdr.mat;
    source_obj_seg(j).hdr.private.mat0=source_obj_seg(j).hdr.mat;
end

big_obj=join_obj(target_object_seg,source_obj_seg);
big_obj=spm_reslice_at(big_obj,wrtflg);%inputs need to be target,source,other
new_source=big_obj(2:end);


end