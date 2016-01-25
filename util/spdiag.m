function mat=spdiag(vec)
sv=size(vec, 1);
mat=sparse([1:sv]', [1:sv]', vec, sv, sv);