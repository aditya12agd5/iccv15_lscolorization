function [left, right, spflag] = spliterrorongmean(leaf, egsmat, eggv, fscales, params)

wmaxleaf = params.maxsplitloc;
minnode = params.minnode;

% for the moment, one random feature
min_split_value=0;
min_split_cost = 1e7;
min_split_projection = [];

grs_record = zeros(params.ntos, wmaxleaf);

if leaf.nindices>2.5*minnode
  minnode_bigger = floor(1.1*minnode);
    for splitc=1:params.ntos
        % split on a random projection of features.  Exp to push it towards chosing a small number.
        projection = randn(params.nfeats_split,1);
        projection = projection.*exp(projection*1.5);
        projection = projection/sum(abs(projection));
        projection = projection.*fscales;
        if(sum(egsmat(leaf.inda,:)*projection == mode(egsmat(leaf.inda,:)*projection)) > .05*numel(leaf.inda));
          continue;
        end

        stepsize = floor((leaf.nindices-2*minnode)/wmaxleaf);
        if(stepsize < 1)
          stepsize = 1;
        end

        leaf_data_indices=leaf.inda;
        funs=egsmat(leaf_data_indices, :)*projection;
        [bvals, bsix]=sort(funs); % this is ascending in the split variable

        all_sorted_indices = leaf_data_indices(bsix);

        l_val = sum(eggv(all_sorted_indices(1:minnode_bigger),:));
        r_val = sum(eggv(all_sorted_indices(minnode_bigger+1:end),:));
        lcount = minnode_bigger;
        rcount = size(all_sorted_indices,1)-minnode_bigger;

        count = 1;
        prev_index = minnode_bigger;
        for index = minnode_bigger:stepsize:size(all_sorted_indices,1)-(minnode_bigger)
          if(prev_index+1 ~= index)
            change = sum(eggv(all_sorted_indices(prev_index+1:index),:));
          else
            change = eggv(all_sorted_indices(prev_index+1:index),:);
          end

          lcount = lcount + stepsize;
          rcount = rcount - stepsize;
          l_val = l_val + change;
          r_val = r_val - change;
          left_mean = l_val / lcount;
          right_mean = r_val / rcount;

          % TEST CODE
          %left_mean_ = mean(eggv(all_sorted_indices(1:index),:));
          %right_mean_ = mean(eggv(all_sorted_indices(index+1:end),:));
          %assert(all(abs(left_mean-left_mean_) < 10^-5));
          %assert(all(abs(right_mean-right_mean_) < 10^-5));

          grs = -(sum(left_mean.^2) * lcount + sum(right_mean.^2) * rcount);
          grs_record(splitc, count) = grs;

          count = count + 1;
          prev_index = index;

          % If the current split beats the previous best.
          if (grs < min_split_cost)
            min_split_cost = grs;
            min_split_value = bvals(index);
            min_split_projection = projection;
          end
        end
    end
end
if isempty(min_split_projection)
    spflag=0;
    left=0;
    right=0;
else
    spflag=1;
    wv=min_split_value+eps;
    nind = egsmat(leaf.inda, :) * min_split_projection < wv;
    lia=leaf.inda(nind==1, :);
    ria=leaf.inda(nind==0, :);

    %disp(['l:' num2str(size(lia,1)) 'r:' num2str(size(ria,1))]);

    if size(lia, 1) < minnode
      fprintf(1, ['l*' num2str(size(lia,1)) 'r*' num2str(size(ria,1)) '-' num2str(minnode) ' ']);
    end
    if size(ria, 1) < minnode
      fprintf(1, ['r*' num2str(size(ria,1)) 'l*' num2str(size(lia,1)) '-' num2str(minnode) ' ']);
    end

    left = node([], [], lia,...
        [leaf.cumsplitprojection, min_split_projection],...
        [leaf.cumsplitval, wv], [leaf.cumsplitsign, -1]);
    %left is always v<val, right is v>val
    % so left is v-val < 0, right is v-val>0
    right = node([], [], ria, ...
        [leaf.cumsplitprojection, min_split_projection],...
        [leaf.cumsplitval, wv], [leaf.cumsplitsign, +1]);
end

