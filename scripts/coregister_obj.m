function [new_source,M]=coregister_obj(target_object_seg,source_obj_seg,flags)
%target object seg should be a dwi object with one element
% operates on "objects" (structures)

%flags: no0: if 1 (0 default), chooses to use a modified spm_coreg that doesn''t consider voxels in first bin of 256 (closest to 0) in joint histogram similarity metric, ie voxels with very low intensity in either source or target will not contribute towards similarity metric (NMI) and therefore registration.

%     : animal: if 1 (0 default), x10 scaled down parameters will be used for spm_coreg. This is very important if your voxel sizes are small (e.g. closer to .1x.1.x.1 than 1.5x1.5x1.5).

%     : apply: if 2 (default), affine registration will be estimated and
%     image will be resliced thusly. If 1, reg will be estimated and the
%     orientation matrix will be updated to reflect the change of origin.
%     If 0, "new_source" will be equal to source_obj_seg and M will be the
%     estimated rigid body reg. 

%     If -1, source_obj_seg will be resliced but
%     not moved to match target_object_seg (Note that this does not mean
%     that target_object_seg won't be moved at all. it will be moved by
%     target.hdr.mat\source.hdr.mat -- aligning the coordinate systems. I
%     have a fix for this that should be an option of this script. I think
%     you can reslice without moving by using align_tool no CoM no reslice
%     and then reslicing (apply=-1) with this script. I'll make sure and
%     maybe add this in the future.

%     : stringent: if 1 (0 default), changes a few options to try and get a
%     better registration with a time cost. No guarantee that registration
%     will be any better

% In the future, normalize how things are output -- currently, a side
% effect of reslicing is that new_source takes on empty fields to match
% target_object_seg. And this doesn't happen if reslicing is not
% performed. This should not practically affect function, but it
% could trip up user-written code.

if length(target_object_seg)>1; error('unintended usage');end

if ~exist('flags','var') || ~isfield(flags,'apply')
    flags.apply=2;
end

if ~exist('flags','var') || ~isfield(flags,'stringent')
    flags.stringent=0;
end

% coregistration and reslicing parameters
estflg.cost_fun = 'nmi';
if flags.stringent~=1
    estflg.sep      = [4 2];
    estflg.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
else %we want more stringent
    estflg.sep      = [2 1]; %sample more 
    estflg.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.0005 0.0005]; %require closer final match
end
    
estflg.fwhm     = [7 7];
wrtflg        = spm_get_defaults('realign.write');
wrtflg.interp   = 1;
wrtflg.which    = [1 0];

if ~isequal(flags.apply,-1) % if the user doesn't want to completely skip coreg and get to reslicing
    
    if isfield(flags,'animal') && isequal(1,flags.animal)
        estflg.sep=estflg.sep*0.1;
        estflg.fwhm=estflg.fwhm*0.1;
    end
    
    
    target_object_seg.hdr.img=target_object_seg.img;
    for i=1:length(source_obj_seg)
        source_obj_seg(i).hdr.img=source_obj_seg(i).img;
    end
    
    
    if ~isfield(flags,'no0') || isequal(flags.no0,0)
        x  = spm_coreg_at(target_object_seg.hdr, source_obj_seg(1).hdr, estflg);
    else
        x  = spm_coreg_at_no0(target_object_seg.hdr, source_obj_seg(1).hdr, estflg);
    end
    
    %unset .hdr.img for clarity
    target_object_seg.hdr=rmfield(target_object_seg.hdr,'img');
    for i=1:length(source_obj_seg)
        source_obj_seg(i).hdr=rmfield(source_obj_seg(i).hdr,'img');
    end
    
    M  = inv(spm_matrix(x));

    if flags.apply>0
        
        MM = zeros(4, 4, length(source_obj_seg));
        
        for j=1:length(source_obj_seg)
            MM(:,:,j)=source_obj_seg(j).hdr.mat;
        end
        for j = 1:length(source_obj_seg)
            source_obj_seg(j).hdr.mat=M*MM(:,:,j);
            source_obj_seg(j).hdr.private.mat=source_obj_seg(j).hdr.mat;
            source_obj_seg(j).hdr.private.mat0=source_obj_seg(j).hdr.mat;
        end
        
        if isequal(flags.apply,1) %if you only want hdr to change, bail out here
            new_source=source_obj_seg;
            return;
        end
    else %if flags.apply is 0, i.e. you just wanted to see what M was, leave source unchanged and bail out here
        new_source=source_obj_seg;
        return
    end
end % endif flags.apply~=-1

%implicit if: other options have bailed at this point
%if ~isfield(flags,'apply') || isequal(flags.apply,2)

%reslice
[big_obj,~,new_src_fields]=join_obj(target_object_seg,source_obj_seg); %this is how spm wants em
big_obj=spm_reslice_at(big_obj,wrtflg);%inputs need to be target,source,other
new_source=big_obj(2:end);
new_source=rmfield(new_source,new_src_fields);

%end implicit "if"
end
