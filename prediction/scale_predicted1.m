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