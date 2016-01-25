function mat=makeprojmat(nums, numt)
% we need to take every second pixel
nipix=nums*numt;
onums=floor((nums+1)/2);
onumt=floor((numt+1)/2);
nopix=onums*onumt;
ivals=[1:nopix]';
jvals=zeros(nopix, 1);
vvals=ones(nopix, 1);
wi=1;
for wy=1:onumt
    for wx=1:onums
        iwx=(wx-1)*2+1;
        iwy=(wy-1)*2+1;
        jvals(wi)=(iwy-1)*nums+iwx;
        wi=wi+1;
    end
end
mat=sparse(ivals, jvals, vvals, nopix, nipix);
    