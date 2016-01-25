function ptr=ijtoptr(i, j, ix)
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
ptr=(j-1)*ix+i;