classdef node < handle

  properties
    % Tree
    leaf;
    left;
    right;

    % Splitting
    cumsplitprojection;
    cumsplitval;
    cumsplitsign;

    % Learning
    inda;
    weights;

    % Model
    amat;
    bvmat;
    alpha;
  end

  methods
    function obj = node(left, right, inda, cumsplitprojection, cumsplitval, cumsplitsign)
      obj.left = left;
      obj.right = right;
      obj.inda = inda;
      obj.cumsplitval = cumsplitval;
      obj.cumsplitsign = cumsplitsign;
      obj.cumsplitprojection = cumsplitprojection;
      obj.leaf = false;
    end

    function indexsize = nindices(obj)
      indexsize = numel(obj.inda);
    end

    function flat = flatten(obj)
      if(obj.leaf)
        flat = obj;
      else
        flat = [obj.left.flatten(), obj.right.flatten()];
      end
    end

    function weights = getWeights(obj, egsmat, fscales)
      weights = ones(size(egsmat,1), 1);
      smoothing_func = @(x) [x.^0 x x.^2 x.^3] * [.5;.75;0;-.25];
      for i = 1:size(obj.cumsplitprojection, 2)
        uv = egsmat*(obj.cumsplitprojection(:,i)) - obj.cumsplitval(i);
        suv = (uv * obj.cumsplitsign(i));
        selection = abs(suv) > 1;
        suv(selection) = sign(suv(selection));
        weights = weights .* smoothing_func(suv);
      end
     end

  end
end
