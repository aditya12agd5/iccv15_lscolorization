function [axmat, aymat] = setupfirstderivmat(nums, numt, dx, dy)
% sparse matrix representations of derivatives

nxcs=nums*numt;
xiv=zeros(2*nxcs, 1);
xjv=zeros(2*nxcs, 1);
xvv=zeros(2*nxcs, 1);
xp=1;
yiv=zeros(2*nxcs, 1);
yjv=zeros(2*nxcs, 1);
yvv=zeros(2*nxcs, 1);
yp=1;

ptr=1;
for i=1:nums
  for j=1:numt
    ptr=(j-1)*nums+i;
    if i>1&i<nums
        xiv(xp)=ptr;
        xjv(xp)=(j-1)*nums+i;
        xvv(xp)=-1/dx;
        xp=xp+1;
        xiv(xp)=ptr;
        xjv(xp)=(j-1)*nums+i+1;
        xvv(xp)=1/dx;
        xp=xp+1;
%       axmat(ptr, (j-1)*nums+i)=-1/(dx);
%       axmat(ptr, (j-1)*nums+i+1)=1/(dx);
%       axxmat(ptr, (j-1)*nums+i-1)=1/dx^2;
%       axxmat(ptr, (j-1)*nums+i)=-2/dx^2;
%       axxmat(ptr, (j-1)*nums+i+1)=1/dx^2;
    elseif i>1
        xiv(xp)=ptr;
        xjv(xp)=(j-1)*nums+i-1;
        xvv(xp)=-1/dx;
        xp=xp+1;
        xiv(xp)=ptr;
        xjv(xp)=(j-1)*nums+i;
        xvv(xp)=1/dx;
        xp=xp+1;
%       axmat(ptr, (j-1)*nums+i-1)=-1/dx;
%       axmat(ptr, (j-1)*nums+i)=1/dx;
%       axxmat(ptr, (j-1)*nums+i-2)=1/dx^2;
%       axxmat(ptr, (j-1)*nums+i-1)=-2/dx^2;
%       axxmat(ptr, (j-1)*nums+i)=1/dx^2;
    elseif i<nums
        xiv(xp)=ptr;
        xjv(xp)=(j-1)*nums+i;
        xvv(xp)=-1/dx;
        xp=xp+1;
        xiv(xp)=ptr;
        xjv(xp)=(j-1)*nums+i+1;
        xvv(xp)=1/dx;
        xp=xp+1;
%       axmat(ptr, (j-1)*nums+i)=-1/dx;
%       axmat(ptr, (j-1)*nums+i+1)=1/dx;
%       axxmat(ptr, (j-1)*nums+i)=1/dx^2;
%       axxmat(ptr, (j-1)*nums+i+1)=-2/dx^2;
%       axxmat(ptr, (j-1)*nums+i+2)=1/dx^2;

    end
    if j>1&j<numt
        yiv(yp)=ptr;
        yjv(yp)=(j-1)*nums+i;
        yvv(yp)=-1/dy;
        yp=yp+1;
        yiv(yp)=ptr;
        yjv(yp)=(j)*nums+i;
        yvv(yp)=1/dy;
        yp=yp+1;
%       aymat(ptr, (j-1)*nums+i)=-1/(dy);
%       aymat(ptr, (j)*nums+i)=1/(dy);
%       ayymat(ptr, (j-2)*nums+i)=1/dy^2;
%       ayymat(ptr, (j-1)*nums+i)=-2/dy^2;
%       ayymat(ptr, (j)*nums+i)=1/dy^2;
    elseif j>1
        yiv(yp)=ptr;
        yjv(yp)=(j-2)*nums+i;
        yvv(yp)=-1/dy;
        yp=yp+1;
        yiv(yp)=ptr;
        yjv(yp)=(j-1)*nums+i;
        yvv(yp)=1/dy;
        yp=yp+1;
%       aymat(ptr, (j-2)*nums+i)=-1/dy;
%       aymat(ptr, (j-1)*nums+i)=1/dy;
%       ayymat(ptr, (j-3)*nums+i)=1/dy^2;
%       ayymat(ptr, (j-2)*nums+i)=-2/dy^2;
%       ayymat(ptr, (j-1)*nums+i)=1/dy^2;
    elseif j<numt
        yiv(yp)=ptr;
        yjv(yp)=(j-1)*nums+i;
        yvv(yp)=-1/dy;
        yp=yp+1;
        yiv(yp)=ptr;
        yjv(yp)=(j)*nums+i;
        yvv(yp)=1/dy;
        yp=yp+1;
%       aymat(ptr, (j-1)*nums+i)=-1/dy;
%       aymat(ptr, (j)*nums+i)=1/dy;
%       ayymat(ptr, (j-1)*nums+i)=1/dy^2;
%       ayymat(ptr, (j)*nums+i)=-2/dy^2;
%       ayymat(ptr, (j+1)*nums+i)=1/dy^2;
    end
  end
end
axmat=sparse(xiv(1:xp-1), xjv(1:xp-1), xvv(1:xp-1), nxcs, nxcs);
aymat=sparse(yiv(1:yp-1), yjv(1:yp-1), yvv(1:yp-1), nxcs, nxcs);
