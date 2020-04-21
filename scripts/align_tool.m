function [new_source,shift,trg_CoM]=align_tool(target_object_seg,source_object_seg,flags)
%BIG NOTE: RIGHT NOW, THIS IS ONLY INTENDED TO BE USED ON images with the
%same dimensions (eg two instances of the same sequence)

%right now, flags.flip does nothing, images are flipped anyway. To use it
%another way probably wouldn't work in most cases.

%target object seg should be a dwi object with one element. 
if length(target_object_seg)>1; error('unintended usage, target object should be length 1');end

if ~exist('flags','var')
    flags=[];
end
if ~isfield(flags,'com')
    flags.com=0;
end
if ~isfield(flags,'reslice')
    flags.reslice=1;
end
if ~isfield(flags,'flip')
    flags.flip=0;
end

%ensure correct usage
assert(all(target_object_seg.hdr.dim == source_object_seg(1).hdr.dim))


%% do flips if user wants to

% change image to match orientation convention of target
% i.e. only do 180 flips along any axis

% this code was just added from shih_coreg_thursday (...) 

mat_flip=eye(4);
if target_object_seg(1).hdr.mat(1)*source_object_seg(1).hdr.mat(1)<0
    %flip
    for i=1:length(source_object_seg)
        source_object_seg(i).img=flip(source_object_seg(i).img,1);
%         mat_flip(1)=-1;
    end
end
if target_object_seg(1).hdr.mat(6)*source_object_seg(1).hdr.mat(6)<0
    %flip
    for i=1:length(source_object_seg)
        source_object_seg(i).img=flip(source_object_seg(i).img,2);
%         mat_flip(6)=-1;
    end
end
if target_object_seg(1).hdr.mat(11)*source_object_seg(1).hdr.mat(11)<0
    %flip
    for i=1:length(source_object_seg)
        source_object_seg(i).img=flip(source_object_seg(i).img,3);
%         mat_flip(11)=-1;
    end
end

%
% for i=1:length(source_object_seg)
%     source_object_seg(i).hdr.mat=mat_flip*source_object_seg(i).mat; 
% end



%% move to align CoM if user has specified so

shift=eye(4); %to be used if ~use_CoM
if flags.com
    trg_CoM=centerOfMass(target_object_seg.img);
    src_CoM=centerOfMass(source_object_seg(1).img);
    
    %create 4x4 of required shift
    src_shift=trg_CoM-src_CoM;
    shift=spm_matrix(src_shift);
end

trg_space=target_object_seg.hdr.mat;
for j = 1:length(source_object_seg)
    %next line causes (in spm_reslice) source voxels to shift, then source
    %hdr.mat to be set to target
    source_object_seg(j).hdr.mat=shift*trg_space;
    
    %set the others to the same
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