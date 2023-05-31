function [bp,p] = xiugai(X, Y, M, varargin) 
    i=1;
    % Defaults

    targetu = .20;              % proportion contribution of boot procedure to p-value
    doplots = 0;                % make plots
    boottop = 0;                  % bootstrap instead of OLS
    dobootfirstlevel = 0;       % bootstrap first level in multi-level analysis
    dorobust = 0;               % robust IRLS
    verbose = 0;                % verbose output
    vnames = {'X' 'Y' 'M'};     % variable names
    bootsamples = 1000;         % initial bootstrap samples
    whpvals_for_boot = [1 2 3 5];   % indices of p-values, the min of which is used to determine boot samples needed
    % lower p-values require more boot samples
    % for the p-vals to be meaningful.
    domultilev = 1;             % multi-level analysis (if N replications > 1)
    arorder = 0;                % timeseries AR(p) model order at first level
    shiftrange = [];            % min and max to shift timeseries by, in elements (or [] for no shift)
    dolatent = 0;               % latent HRF model
    doCIs = 0;

    dosave = 0;                 % save figures at end

    dosignperm = 0;             % sign permutation (alt. to bootstrap)
    persistent_perms = 0;       % keep same permutation matrix across repeated calls to mediation

    additionalM = [];           % additional mediators, entered in same format as other inputs

    X_2ndlevel = [];
    l2mnames = {};              % Level-2 moderator names
    logistic_Y = 0;
    [paths(i,:), abetas{i}, bbetas{i}, cpbetas{i}, cbetas{i}, sterrs(i,:), intcpt, n(i), residm{i}, residy{i}, residy2{i}, ...
                Vxm(i, :), Vmy(i), Vxy(i)] = ...
                mediation_path_coefficients(X, Y, M, domultilev, dorobust, boottop, logistic_Y);
%     a=abetas{i}(1,1);
%     b=bbetas{i}(1,1);
% n为时间点个数
    bp=paths(i,:);
    ste=sterrs(i,:);   %标准差
    
    t = bp ./ ste;  %检验值
    df = [n-1 n-2 n-2 n-1 n-2];  %自由度
    p = min(1, (2 .* (1 - tcdf(abs(t), df))));  %% tcdf：t分布的累计概率密度函数 计算t检验的p值
%     b_p=p(1,2);
%     cp=cpbetas{i}(1,1);
%     cp_p=p(1,3); 
%     ab = a .* b;
%     ab_p=p(1,5);
end

