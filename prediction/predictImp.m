pa='H:\Preprocessing\HCP360\REST1\result\med\predictModel\liangbiao12.mat';
load(pa);
[n_sub, n_scale] = size(liangbiaoAll);
thresh = 0.05;
M = nan(n_scale, 2);
N_MA = nan(n_sub, n_scale);
N_MS = nan(n_sub, n_scale);

for i = 1:n_scale
    liangbiao = liangbiaoAll(:, i);
    [M(i, 1), N_MA(:, i)] =  scale_predicted1(A, liangbiao, age, sex, thresh);
    [M(i, 2), N_MS(:, i)] =  scale_predicted1(B, liangbiao, age, sex, thresh);
end
