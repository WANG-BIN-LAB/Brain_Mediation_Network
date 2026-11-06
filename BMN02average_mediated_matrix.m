% Average across valid mediators for each ordered source-target pair.
n = size(tensor,1);
assert(n == size(tensor,2) && n == size(tensor,3), ...
    'Expected a source x target x mediator tensor.');
AMM = mean(tensor,3,'omitnan');
