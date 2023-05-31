pa='H:\Preprocessing\HCP360\REST1\result\*.mat';
Listn =dir(pa);    %被试的输入，到病人正常人最外面那个路径
   %所有的cell(i,1)相加求平均得到264x264的矩阵
  %反应每个脑区被中介的能力  之后求矩阵某一节点该列加该行的平均值
AVE_Sum=zeros(360,360);  
 
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
    savefile='H:\Preprocessing\HCP360\REST1\result\BeiMedMAX\';
    save(strcat(savefile,'M_',num2str(n),'_Med'),'AVE'); 
end
