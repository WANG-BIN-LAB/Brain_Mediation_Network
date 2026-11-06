function [R, behav_pred, fold_info] = cpm_loocv_nested(all_mats, all_behav, age, sex, thresh)
% all_mats: nodes x subjects (MA or MS). Encode sex as 0/1.
% All vectors must follow the column order of all_mats.
if nargin < 5
    thresh = 0.05;
end
validateattributes(all_mats, {'numeric'}, {'2d','nonempty','real','finite'});
X = double(all_mats');
n = size(X, 1);
validateattributes(all_behav, {'numeric'}, {'vector','numel',n,'real','finite'});
validateattributes(age, {'numeric'}, {'vector','numel',n,'real','finite'});
validateattributes(sex, {'numeric','logical'}, {'vector','numel',n,'binary'});
validateattributes(thresh, {'numeric'}, {'scalar','real','>',0,'<',1});
assert(n >= 5, 'At least five subjects are required.');

y = double(all_behav(:));
C = [ones(n,1), double(age(:)), double(sex(:))];
behav_pred = nan(n,1);
fold_info.n_pos = zeros(n,1);
fold_info.n_neg = zeros(n,1);

for held_out = 1:n
    train = true(n,1);
    train(held_out) = false;
    C_train = C(train,:);
    C_test = C(held_out,:);
    assert(rank(C_train) == 3, 'Age and sex must be estimable in each training fold.');

    x_coef = C_train \ X(train,:);
    y_coef = C_train \ y(train);
    x_train = X(train,:) - C_train * x_coef;
    x_test = X(held_out,:) - C_test * x_coef;
    y_train = y(train) - C_train * y_coef;

    usable = std(x_train,0,1) > sqrt(eps) .* max(abs(X(train,:)),[],1);
    x_train = x_train(:,usable);
    x_test = x_test(:,usable);
    pos = false(1,size(x_train,2));
    neg = pos;
    if any(usable) && norm(y_train) > sqrt(eps) * norm(y(train))
        [rho,p] = corr(x_train,y_train,'Type','Spearman');
        pos = (p < thresh & rho > 0)';
        neg = (p < thresh & rho < 0)';
    end
    fold_info.n_pos(held_out) = nnz(pos);
    fold_info.n_neg(held_out) = nnz(neg);

    train_scores = [sum(x_train(:,pos),2), sum(x_train(:,neg),2)];
    test_scores = [sum(x_test(:,pos),2), sum(x_test(:,neg),2)];
    center = mean(train_scores,1);
    scale = std(train_scores,0,1);
    active = scale > sqrt(eps) .* max(abs(train_scores),[],1);
    train_z = (train_scores(:,active) - center(active)) ./ scale(active);
    test_z = (test_scores(:,active) - center(active)) ./ scale(active);

    beta = lsqminnorm([ones(nnz(train),1), train_z], y_train);
    % Restore the age/sex baseline to predict on the original scale.
    behav_pred(held_out) = C_test * y_coef + [1, test_z] * beta;
end

fold_info.n_feat = fold_info.n_pos + fold_info.n_neg;
R = corr(behav_pred,y,'Type','Spearman');
end
% Copyright 2015 Xilin Shen and Emily Finn 

% This code is released under the terms of the GNU GPL v2. This code
% is not FDA approved for clinical use; it is provided
% freely for research purposes. If using this in a publication
% please reference this properly as: 

% Finn ES, Shen X, Scheinost D, Rosenberg MD, Huang, Chun MM,
% Papademetris X & Constable RT. (2015). Functional connectome
% fingerprinting: Identifying individuals using patterns of brain
% connectivity. Nature Neuroscience 18, 1664-1671.

% This code provides a framework for implementing functional
% connectivity-based behavioral prediction in a leave-one-subject-out
% cross-validation scheme, as described in Finn, Shen et al 2015 (see above
% for full reference). The first input ('all_mats') is a pre-calculated
% MxMxN matrix containing all individual-subject connectivity matrices,
% where M = number of nodes in the chosen brain atlas and N = number of
% subjects. Each element (i,j,k) in these matrices represents the
% correlation between the BOLD timecourses of nodes i and j in subject k
% during a single fMRI session. The second input ('all_behav') is the
% Nx1 vector of scores for the behavior of interest for all subjects.

% As in the reference paper, the predictive power of the model is assessed
% via correlation between predicted and observed scores across all
% subjects. Note that this assumes normal or near-normal distributions for
% both vectors, and does not assess absolute accuracy of predictions (only
% relative accuracy within the sample). It is recommended to explore
% additional/alternative metrics for assessing predictive power, such as
% prediction error sum of squares or prediction r^2.

function [R, P,behav_pred]=scale_predicted(percontro_BD,trancontro_BD,SAPS)
% clear;
% clc;

% ------------ INPUTS -------------------

% all_mats  = contro_mol_BP;
all_behav = SAPS;

% threshold for feature selection
thresh = 0.05;
[r_mat,p_mat] = corr(percontro_BD',all_behav,'type','Spearman');
 node = find( p_mat < thresh);
 sign=node;
 all_mats1=percontro_BD(sign,:);
 [r_mat,p_mat] = corr(trancontro_BD',all_behav,'type','Spearman');
 node = find( p_mat < thresh);
 sign=node;
 all_mats2=trancontro_BD(sign,:);
 
%  [r_mat,p_mat] = corr(contro_syn_BP',all_behav,'type','Spearman');
%  node = find( p_mat < thresh);
%  sign=node;
%  all_mats3=contro_syn_BP(sign,:);
 all_mats=[all_mats1;
     all_mats2];

no_sub = size(all_mats,2);
% no_node = size(all_mats,1);
% 
% behav_pred_pos = zeros(no_sub,1);
% behav_pred_neg = zeros(no_sub,1);
behav_pred = zeros(no_sub,1);
% sign=node;
% all_mats=all_mats(sign,:);
no_node = size(all_mats,1);
for leftout = 1:no_sub;
    fprintf('\n Leaving out subj # %6.3f',leftout);
    
    % leave out subject from matrices and behavior
    
    train_mats = all_mats;
    train_mats(:,leftout) = [];
    train_vcts = train_mats;
    
    train_behav = all_behav;
    train_behav(leftout) = [];
    
    % correlate all edges with behavior

    [r_mat,p_mat] = corr(train_vcts',train_behav,'type','Spearman');
    
%     r_mat = reshape(r_mat,no_node,no_node);
%     p_mat = reshape(p_mat,no_node,no_node);
    
    % set threshold and define masks
    
    pos_mask = zeros(no_node,1);
    neg_mask = zeros(no_node,1);
    
    pos_edges = find(r_mat > 0 & p_mat < thresh);
    neg_edges = find(r_mat < 0 & p_mat < thresh); 
    
    pos_mask(pos_edges) = 1;
    neg_mask(neg_edges) = 1;
    
    % get sum of all edges in TRAIN subs (divide by 2 to control for the
    % fact that matrices are symmetric)
    
    train_sumpos = zeros(no_sub-1,1);
    train_sumneg = zeros(no_sub-1,1);
    
    for ss = 1:size(train_sumpos);
        train_sumpos(ss) = sum(sum(train_mats(:,ss).*pos_mask))/2;
        train_sumneg(ss) = sum(sum(train_mats(:,ss).*neg_mask))/2;
    end
    
    % build model on TRAIN subs

%     fit_pos = polyfit(train_sumpos, train_behav,1);
%     fit_neg = polyfit(train_sumneg, train_behav,1);
    b=regress(train_behav,[train_sumpos,train_sumneg,ones(no_sub-1,1)]);
    % run model on TEST sub
     
    test_mat = all_mats(:,leftout);
    test_sumpos = sum(sum(test_mat.*pos_mask))/2;
    test_sumneg = sum(sum(test_mat.*neg_mask))/2;
   
%     behav_pred_pos(leftout) = fit_pos(1)*test_sumpos + fit_pos(2);
%     behav_pred_neg(leftout) = fit_neg(1)*test_sumneg + fit_neg(2);
    behav_pred(leftout) =b(1)*test_sumpos + b(2)*test_sumneg+b(3);
    
end

% compare predicted and observed scores

% [R_pos, P_pos] = corr(behav_pred_pos,all_behav,'type','Spearman')
% [R_neg, P_neg] = corr(behav_pred_neg,all_behav,'type','Spearman')
[R, P] = corr(behav_pred,all_behav,'type','Spearman');

% figure(1); plot(behav_pred_pos,all_behav,'r.'); lsline
% figure(2); plot(behav_pred_neg,all_behav,'b.'); lsline
figure; plot(behav_pred,all_behav,'r.'); lsline
