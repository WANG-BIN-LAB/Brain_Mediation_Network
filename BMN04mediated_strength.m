pa='H:\Preprocessing\HCP360\REST1\result\*.mat';
Listn =dir(pa);    %被试的输入，到病人正常人最外面那个路径
   %所有的cell(i,1)相加求平均得到246x246的矩阵
  %反应每个脑区被中介的能力  之后求矩阵某一节点该列加该行的平均值
AVE_beimed=zeros(500,360);
h=1;
for n=1:500
    FileName = ['B_sub_' num2str(n) '.mat'];
    B=importdata(FileName);
    Q=zeros(360,360); 
    for k=360
       B{k,1}(360,360)=0;
    end
    for k=1:360     
        Q=Q+B{k,1};
    end
    AVE=Q/360;
    for i=1:360
        a=sum(AVE(i,:),2);  %对某一行按列求和
        b=sum(AVE(:,i),1);  %对某一列按行求和
        AVE_beimed(h,i)=(a+b)/720;
    end
    h=h+1;
end
savefile='H:\Preprocessing\HCP360\REST1\result\beimed\all\';
save(strcat(savefile,'beimed_sub2.mat'),'AVE_beimed');