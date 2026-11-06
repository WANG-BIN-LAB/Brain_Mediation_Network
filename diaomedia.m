function [tensor, info] = diaomedia(roi, epsilonT, epsilonI, epsilonD, icTimeCourses, seed)
validateattributes(roi, {'numeric'}, {'2d','real','finite','nonempty'});
validateattributes(epsilonT, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(epsilonI, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(epsilonD, {'numeric'}, {'scalar','real','finite','nonnegative'});
roi = double(roi);
[nTime,nParcels] = size(roi);
assert(nTime > 23 && nParcels >= 20, 'Insufficient time points or parcels for 20 ICs.');

if isempty(icTimeCourses)
    icTimeCourses = bmn_fastica(roi,20,seed);
end
validateattributes(icTimeCourses, {'numeric'}, ...
    {'size',[nTime,20],'real','finite'});
covariates = [ones(nTime,1), double(icTimeCourses)];
assert(rank(covariates) == 21, 'IC time courses and intercept must have full rank.');

% Residualization gives the same path slopes as including the same IC
% covariates and intercept in each of the three OLS equations.
[Q,~] = qr(covariates,0);
residualROI = roi - Q * (Q' * roi);
gram = residualROI' * residualROI;
ss = diag(gram);
validNode = ss > eps .* sum((roi - mean(roi,1)).^2,1)';
total = gram ./ ss;
diagonal = eye(nParcels) ~= 0;

tensor = nan(nParcels,nParcels,nParcels);
counts = zeros(nParcels,5);

for m = 1:nParcels
    eligible = ~diagonal;
    eligible(m,:) = false;
    eligible(:,m) = false;
    possible = eligible & validNode & validNode' & validNode(m);

    alpha = gram(:,m) ./ ss;
    mediatorSS = ss(m) - alpha .* gram(:,m);
    estimable = possible & (mediatorSS > sqrt(eps) * ss(m));
    beta = (gram(m,:) - alpha .* gram) ./ mediatorSS;
    indirect = alpha .* beta;
    direct = total - indirect;
    estimable = estimable & isfinite(indirect) & isfinite(total) & isfinite(direct);

    nearZero = estimable & abs(total) <= epsilonT;
    stable = estimable & ~nearZero;
    noMediation = stable & abs(indirect) <= epsilonI;
    suppression = stable & ~noMediation & ...
        ((indirect > 0 & direct < 0) | (indirect < 0 & direct > 0));
    complete = stable & ~noMediation & ~suppression & abs(direct) <= epsilonD;
    partial = stable & ~noMediation & ~suppression & ~complete;

    values = zeros(nParcels);
    values(complete) = 1;
    values(partial) = indirect(partial) ./ total(partial);
    included = noMediation | complete | partial;
    counts(m,:) = [nnz(eligible), nnz(included), nnz(nearZero), ...
        nnz(suppression), nnz(eligible & ~estimable)];
    values(~included) = NaN;
    tensor(:,:,m) = values;
end
info.tripletCounts = array2table(counts, 'VariableNames', ...
    {'Eligible','Included','NearZeroTotal','Suppression','NonEstimable'});
info.epsilon = struct('total',epsilonT,'indirect',epsilonI,'direct',epsilonD);
info.IC = icTimeCourses;
end
