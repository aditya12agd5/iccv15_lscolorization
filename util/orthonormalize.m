function [A_m] = orthonormalize(A_m, proj_dimensions, nrecon)

% this should be a gram-schmidt orthonormalization of the A
% portion of the A,B,b vector.  It's pretty finicky, but it
% does work.
for i=1:proj_dimensions
    for j=1:i-1
        A_m(:, (i-1)*nrecon+1:i*nrecon)=A_m(:, (i-1)*nrecon+1:i*nrecon)...
            -(sum(A_m(:, (j-1)*nrecon+1:j*nrecon).*A_m(:, (i-1)*nrecon+1:i*nrecon), 2)*ones(1, nrecon)).*...
            A_m(:, (j-1)*nrecon+1:j*nrecon);
    end
    wvs=sum(A_m(:, (i-1)*nrecon+1:i*nrecon).*A_m(:, (i-1)*nrecon+1:i*nrecon), 2).^(1/2);
    if(wvs ~= 0)
      A_m(:, (i-1)*nrecon+1:i*nrecon)=((1./wvs)*ones(1, nrecon)).*A_m(:, (i-1)*nrecon+1:i*nrecon);
    end
end

end
