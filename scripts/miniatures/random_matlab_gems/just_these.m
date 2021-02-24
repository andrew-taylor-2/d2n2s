function filtered_obj=just_these(input_img_object,indices_to_keep,binarize_bool)

%.img fields are assumed to hold a 3d image

if ~exist('binarize_bool','var')
    binarize_bool=0;
end

assert(isvector(indices_to_keep))

for i=1:length(input_img_object)
    good_voxels_4d=abs(input_img_object(i).img-reshape(indices_to_keep,[1 1 1 length(indices_to_keep)]))<0.3; %this implements thresholding if it's a probabilistic region, but this function wasn't intended for probabilistic regions
    good_voxels=logical(sum(good_voxels_4d,4));
    filtered_obj(i)=input_img_object(i);
    filtered_obj(i).img(~good_voxels)=0;
    
    if binarize_bool
        filtered_obj(i).img(good_voxels)=1;
    end
end