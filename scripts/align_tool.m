function [new_source,shift,trg_CoM]=align_tool(target_object_seg,source_object_seg,flags)
%BIG NOTE: RIGHT NOW, THIS IS ONLY INTENDED TO BE USED ON images with the
%same dimensions (eg two instances of the same sequence)

%target object seg should be a dwi object with one element. 
if length(target_object_seg)>1; error('unintended usage');end

if ~exist('flags','var')
    flags=[];
end
if ~isfield(flags,'com')
    flags.com=1;
end
if ~isfield(flags,'reslice')
    flags.reslice=1;
end


%ensure correct usage
% assert(all(target_object_seg.hdr.dim == source_object_seg(1).hdr.dim))


%% move to align CoM if user has specified so

shift=eye(4); %to be used if ~use_CoM
if flags.com
    trg_CoM=centerOfMass(target_object_seg.img);
    src_CoM=centerOfMass(source_object_seg(1).img);
    
    %create 4x4 of required shift
    src_shift=trg_CoM-src_CoM;
    shift=spm_matrix(src_shift);
end

trg_space_params=spm_imatrix(target_object_seg.hdr.mat);



for j = 1:length(source_object_seg)
    %next line causes (in spm_reslice) source voxels to shift, then source
    %hdr.mat to be take on target's origin (but keep its own voxel size)
    
    src_space_params=spm_imatrix(source_object_seg(j).hdr.mat);
    sign_change=arrayfun(@(ii) sign(src_space_params(ii))*sign(trg_space_params(ii)),7:9);
    src_with_trg_origin=spm_matrix([trg_space_params(1:6) src_space_params(7:9).*sign_change trg_space_params(10:12)]);
    
    %set all records of orientation to this just in case
    source_object_seg(j).hdr.mat=shift*src_with_trg_origin;
    source_object_seg(j).hdr.private.mat=source_object_seg(j).hdr.mat;
    source_object_seg(j).hdr.private.mat0=source_object_seg(j).hdr.mat;
end

big_obj=join_obj(target_object_seg,source_object_seg);
if flags.reslice % you don't need to reslice if you're not moving anything
    wrtflg        = spm_get_defaults('realign.write');
    wrtflg.interp   = 1;
    wrtflg.which    = [1 0];

    big_obj=spm_reslice_at(big_obj,wrtflg);%inputs need to be target,source,other
end
new_source=big_obj(2:end);


end