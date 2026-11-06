% For each endpoint, average its available incoming and outgoing AMM entries.
assert(ismatrix(AMM) && size(AMM,1) == size(AMM,2), ...
    'Expected a square AMM matrix.');
MS = mean([AMM, AMM'],2,'omitnan');
