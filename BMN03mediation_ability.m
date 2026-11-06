% For each mediator, average across valid ordered source-target pairs.
n = size(tensor,1);
assert(n == size(tensor,2) && n == size(tensor,3), ...
    'Expected a source x target x mediator tensor.');
MA = mean(reshape(tensor,[],n),1,'omitnan')';
