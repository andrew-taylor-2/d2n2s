function fn = fnify2(dir_object,bhvr_bin)
% bhvr_bin = 1: error when dir_object is empty
% bhvr_bin = 0: fn is empty string when dir_object is empty
if ~exist('bhvr_bin','var')
    bhvr_bin=0;
end
if isempty(dir_object) && bhvr_bin
    error('your dir call was empty, so I can''t make a filename from it')
elseif isempty(dir_object) && ~bhvr_bin
    fn='';
    return
elseif numel(dir_object)>1
    warning(['numel of object from dir call is greater than 1.\n -- it''s ' num2str(numel(dir_object)) ' -- using the first.']);
end

fn=[dir_object(1).folder filesep dir_object(1).name];
end

