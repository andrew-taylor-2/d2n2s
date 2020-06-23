function [out_obj,flip_mat]=flip_mat_obj(in_obj,dim)

flip_mat=eye(4);
if dim==1, ind=1; elseif dim==2, ind=6; elseif dim==3, ind=11; else, error;end
flip_mat(ind)=-1;

in_obj_flipped=in_obj;
in_obj_flipped.hdr.mat=flip_mat*in_obj_flipped.hdr.mat;

out_obj=coregister_obj(in_obj,in_obj_flipped,make_flags('coregister','apply',-1));
end