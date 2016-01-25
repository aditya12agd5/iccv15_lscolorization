function [i, j]=ptrtoij(ptr, ix)
j=floor((ptr-1e-9)/ix)+1;
i=ptr-(j-1)*ix;


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