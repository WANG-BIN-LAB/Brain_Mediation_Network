function [icTimeCourses, iterations] = bmn_fastica(roi, nComponents, seed)
% Symmetric FastICA with tanh nonlinearity.
% Input: time points x ROI. Output: time points x independent components.
if nargin < 2
    nComponents = 20;
end
if nargin < 3
    seed = 1;
end
validateattributes(roi, {'numeric'}, {'2d','real','finite','nonempty'});
validateattributes(nComponents, {'numeric'}, {'scalar','integer','positive'});
validateattributes(seed, {'numeric'}, {'scalar','integer','>=',0,'<=',2^32-1});
X = double(roi);
nTime = size(X,1);
assert(nComponents <= min(nTime-1,size(X,2)), ...
    'Too few time points or parcels for the requested IC count.');

X = X - mean(X,1);
parcelSD = std(X,0,1);
assert(all(parcelSD > 0), 'ICA input contains a constant parcel.');
X = X ./ parcelSD;

% PCA whitening: each retained direction has unit mean squared amplitude.
[U,S,~] = svd(X,'econ');
singularValues = diag(S);
rankTolerance = max(size(X)) * eps(singularValues(1));
assert(singularValues(nComponents) > rankTolerance, ...
    'The standardized ROI matrix has insufficient rank.');
Z = sqrt(nTime) * U(:,1:nComponents);

stream = RandStream('mt19937ar','Seed',seed);
[W,~] = qr(randn(stream,nComponents,nComponents),0);
tolerance = 1e-6;
maxIterations = 1000;
for iterations = 1:maxIterations
    projected = Z * W;
    g = tanh(projected);
    derivativeMean = mean(1 - g.^2,1);
    update = (Z' * g) / nTime - W .* derivativeMean;

    % The orthogonal polar factor gives symmetric decorrelation.
    [left,~,right] = svd(update,'econ');
    nextW = left * right';
    change = max(abs(1 - abs(diag(W' * nextW))));
    W = nextW;
    if change < tolerance
        icTimeCourses = Z * W;
        return
    end
end
error('bmn_fastica:NoConvergence', ...
    'FastICA did not converge in %d iterations. Try a different seed.', maxIterations);
end
