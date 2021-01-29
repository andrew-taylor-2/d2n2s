function filtered_obj=just_these(input_img_object,indices_to_keep)

%.img fields are assumed to hold a 3d image

%this doesnt really need the inputs/outs to be objects. I could just give
%it the .img field and make it more general

assert(isvector(indices_to_keep))

for i=1:length(input_img_object)
    good_voxels_4d=abs(input_img_object(i).img-reshape(indices_to_keep,[1 1 1 length(indices_to_keep)]))<0.3; %this implements thresholding if it's a probabilistic region, but this function wasn't intended for probabilistic regions
    good_voxels=logical(sum(good_voxels_4d,4));
    filtered_obj(i)=input_img_object(i);
    filtered_obj(i).img(~good_voxels)=0;
end