function [out_obj,flip_mat]=flip_mat_obj(in_obj,dim)
% good function

flip_mat=eye(4); %start with identity matrix
if dim==1, ind=1; elseif dim==2, ind=6; elseif dim==3, ind=11; else, error;end %change sign of one of the diagonal elements
flip_mat(ind)=-1;

in_obj_flipped=in_obj;
in_obj_flipped.hdr.mat=flip_mat*in_obj_flipped.hdr.mat; %use what we just made to flip the orientation matrix

out_obj=coregister_obj(in_obj,in_obj_flipped,make_flags('coregister','apply',-1)); %this command essentially does spm_reslice -- the final image will end up with the orient. matrix of the original image (not flipped). Because the "source" and "target" have orientations a 180 degree rotation apart, spm_reslice will flip "source" 180 degrees in the dimension you chose
end