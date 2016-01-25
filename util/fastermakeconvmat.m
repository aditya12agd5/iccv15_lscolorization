function wmat=makeconvmat(kernel, kx, ky, ix, iy)
% recall
% reshape([1:12], 3, 4)
% 
% ans =
% 
%      1     4     7    10
%      2     5     8    11
%      3     6     9    12
%
%  this means that ind=(j-1)*ix+i
%
%    and j=floor(ind/ix)+1 and i=ind-(j-1)*ix
%
% r_ij = sum k_{uv} i_{u-i v-j)
%
% 
nonzeros=kx*ky*ix*iy;
ivec=zeros(nonzeros, 1);
jvec=zeros(nonzeros, 1);
vvec=zeros(nonzeros, 1);
wvptr=1;
kwx=floor(kx/2);
kwy=floor(ky/2);
if (~((kx-(2*kwx+1))==0))||(~((ky-(2*kwy+1))==0))
    error('odd kernel sizes only');
end
for rptr=1:(ix*iy)
    % the result; now make this row
    [wi, wj]=ptrtoij(rptr, ix);
    for u=1:kx
        kcenu=u-kx-1;
        imi=wi-kcenu-2;
        if (imi>0)&&(imi<=ix)
            for v=1:ky
                kcenv=v-ky-1;
                imj=wj-kcenv-3;
                if (imj>0)&&(imj<=iy)
                    ivec(wvptr)=rptr;
                    jvec(wvptr)=(imj-1)*ix+imi; % equivalent to ijtoptr(imi, imj, ix);
                    vvec(wvptr)=kernel(u, v);
                    wvptr=wvptr+1;
                end
            end
        end
    end
end
ivec=ivec(1:wvptr-1);
jvec=jvec(1:wvptr-1);
vvec=vvec(1:wvptr-1);
wmat=sparse(ivec, jvec, vvec, ix*iy, ix*iy);
