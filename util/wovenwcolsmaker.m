function [wcols, patches, weave] = wovenwcolsmaker(nums, numt, params)

% TODO(jrock): It turns out this doesn't work as I expected.  This constructs
% patches in the right order and then weaves the patches as expected.  It
% means in order to incorporate a new non-patch channel I would have to add a
% pass through portion to the patches matrix and a channel-by-channel weave to
% the current patch based weave.

hbw = params.half_patch_width;
nchannels = params.nchannels;

winsize=((2*hbw+1)^2);

npix=nums*numt;
tstep=1;

ntiles=nums*numt;

is=zeros(winsize*ntiles, 1);
js=zeros(winsize*ntiles, 1);

ptr=1;
tileptr=1;

%[ll, kk] = meshgrid(1:2*hbw+1, 1:2*hbw+1);
%ll = ll(:);
%kk = kk(:);
%mm = 1:numel(kk);

for j=1:tstep:numt
  for i=1:tstep:nums
    tcolptr=1;
    for l=1:2*hbw+1
      for k=1:2*hbw+1
        % same as sub2ind(.,.,[nums, numt]) but much faster.
        i_ = (i-hbw-1+k);
        j_ = (j-hbw-1+l-1);

        if(~(i_<=0 || i_ > nums || j_ < 0 || j_ >= numt))
          f_val = i_ + j_ * nums;
          is(winsize*(tileptr-1)+tcolptr)=f_val;
          js(winsize*(tileptr-1)+tcolptr)=winsize*(tileptr-1)+tcolptr;
        end
        % we also want the indices of the relevant pixels
        % in mat
        tcolptr=tcolptr+1;
      end
    end

    % This isn't actually faster that the for loops
    %fvaltot = (i-hbw-1+kk) + (j-hbw-1+ll-1) * nums;
    %is(winsize*(tileptr-1)+1:winsize*(tileptr-1)+numel(fvaltot)) = fvaltot;
    %js(winsize*(tileptr-1)+1:winsize*(tileptr-1)+numel(fvaltot)) = winsize*(tileptr-1)+mm;

    %atiles(ptr, :)=reshape(asi(i-hbw:i+hbw, j-hbw:j+hbw),1, (2*hbw+1)*(2*hbw+1));
    %stiles(ptr, :)=reshape(ssi(i-hbw:i+hbw, j-hbw:j+hbw),1, (2*hbw+1)*(2*hbw+1));
    ptr=ptr+1;
    tileptr=tileptr+1;
  end
end
select = is~=0;
is = is(select);
js = js(select);

% this takes a single nums*numt vector, and makes a bunch of winsize x winsize tiles
patches=(sparse(is, js, ones(numel(is),1), npix, winsize*ntiles))';

iw = zeros(nchannels*winsize*ntiles, 1);
jw = zeros(nchannels*winsize*ntiles, 1);
tptr = 1;
[wi, wj, wv]=find(eye(winsize));
for i=1:ntiles
  for j=1:nchannels
    iw(tptr:tptr+winsize-1)=wi+(i-1)*nchannels*winsize + (j-1) * winsize;
    jw(tptr:tptr+winsize-1)=wj+(i-1)*winsize + (j-1) * ntiles * winsize;
    tptr=tptr+winsize;
  end
end

weave=sparse(iw, jw, ones(nchannels*winsize*ntiles, 1), nchannels*winsize*ntiles, nchannels*winsize*ntiles);
wcols=weave*kron(eye(nchannels), patches);

