function spm=spid(n)
spm=sparse([1:n]', [1:n]', ones(n, 1), n, n);