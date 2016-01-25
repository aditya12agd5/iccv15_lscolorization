function [I_bmy_rmg] = converttobmyrmg(I)
  eps = 0.0001;
  gray = sum(I,3)./3;
  I_bmy_rmg= double(zeros(size(I,1), size(I,2), 2));
  I_bmy_rmg(:,:,1) = (I(:,:,3)-(0.5.*(I(:,:,1)+I(:,:,2))))./(gray+eps); 
  I_bmy_rmg(:,:,2) = (I(:,:,1)-I(:,:,2))./(gray+eps);
end
