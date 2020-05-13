function inputimg=binarize_to_value(inputimg,valuee,thresh)
if ~exist('thresh','var')
    thresh=0.01;
end
if ~exist('valuee','var')
    valuee=1;
end
inputimg(inputimg<thresh)=0;
inputimg(inputimg>=thresh)=valuee;
end