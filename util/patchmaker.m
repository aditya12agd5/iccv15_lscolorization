function [patches] = patchmaker(nums, numt, shifts, shiftt, hbw)

winsize=((2*hbw+1)^2);

% In case we don't want to place a patch at every pixel
ntiles = nums*numt;
npix = nums*numt;

is=zeros(winsize*ntiles, 1);
js=zeros(winsize*ntiles, 1);

tileptr=1;
tstep = 1;

for j=1:tstep:numt
  for i=1:tstep:nums
    indptr=1;
    tcolptr=1;
    for l=1:2*hbw+1
      for k=1:2*hbw+1
        i_ = i + shifts - hbw - 1 + k;
        j_ = j + shiftt - hbw - 1 + l - 1;

        if(~(i_ <= 0 || i_ > nums || j_ < 0 || j_ >= numt))
          % same as sub2ind(.,.,[nums, numt]) but much faster.
          f_val = i_ + j_ * nums;
          if(f_val <= 0 || f_val > npix)
            error('fval <=0');
          end
          is(winsize*(tileptr-1)+tcolptr)=f_val;
          js(winsize*(tileptr-1)+tcolptr)=winsize*(tileptr-1)+tcolptr;
        end

        tcolptr=tcolptr+1;
      end
    end
    tileptr=tileptr+1;
  end
end
select = is~=0;
is = is(select);
js = js(select);

% this takes a single nums*numt vector, and makes a bunch of hbw x hbw tiles
patches = (sparse(is, js, ones(numel(is), 1), npix, winsize*ntiles))';
