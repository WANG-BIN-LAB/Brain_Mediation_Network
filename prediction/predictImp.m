pa='H:\Preprocessing\HCP360\REST1\result\med\predictModel\liangbiao12.mat';
load(pa);
M=zeros(12,2);
N=zeros(500,12);
for i = 1:12
    liangbiao=liangbiaoAll(:,i);
    [R, P,behav_pred] = scale_predicted1(B,A,liangbiao);
    N(:,i) = behav_pred;
    M(i,1) = R;
    M(i,2) = P;
end
