function vec=applymat(amat, deps, vec)
vec = amat'* (amat * vec) + deps * vec;
%vec = wcols'*(nwmat'*(nwmat*(wcols*vec)))+deps*vec;
% TODO(jrock): not sure what to make of this.
% (w'm'mw)*v + d*v
