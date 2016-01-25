function spec = createImageSpec(num_indices, func, name)
current_ind = 1;

% setup func
if(nargin < 2)
    func = @no_op;
end
assert(length(func) == length(num_indices) || length(func) == 1);

if(length(func) == 1)
    func_ = func;
    func = {};
    for i = 1:length(num_indices)
        func{i} = func_;
    end
end

% Setup name
if(nargin < 3)
    name = ''
end

assert((iscellstr(name) && length(name) == length(num_indices)) || ischar(name));
if(ischar(name))
    name_ = name;
    name = {};
    for i = 1:length(num_indices)
        if(~strcmp(name_, ''))
            name{i} = [name_ '_' num2str(i)];
        else
            name{i} = '';
        end
    end
end


for i = 1:length(num_indices)
    spec(i) = struct('sample',current_ind:current_ind+num_indices(i)-1,'func',func{i},'name',name{i});
    current_ind = current_ind + num_indices(i);
end

end

function x = no_op(x)
end
