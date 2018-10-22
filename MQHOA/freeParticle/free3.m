% function freeParticle
clear all;         %清除所有变量，包括全局变量
clf;
%hold on;
feature jit off;   %关闭jit加速器
format long;       %控制命令窗口显示方式和位数
global DIM;        %定义全局变量DIM，可以为本源文件中其它函数所公用，有效范围为从定义变量的位置开始到本源文件结束。
DIM=1;             %目标函数的维度
repeat=30;          %重复计算的次数  
sigmaMin=0.000001;  %最小尺度
gap=0.00001;   %稳态差值
groupNum=12;        %高斯采样的数目，k值
minDomain=-5;maxDomain=5; %目标函数的定义域上下界  
gbestV=zeros(1,repeat);      %每次计算的最优函数值
gfe=zeros(1,repeat);         %每次计算的函数进化次数
tot_time = zeros(1,repeat);  %每次计算的时间
%MFE=10000;  %最大迭代次数
badspot=0;
tic;           %计时开始  
funcV=1./zeros(1,groupNum);        %存储k个采样点的函数值（向量）
samplePos=zeros(DIM,groupNum);  %当前k个采样点坐标（矩阵）
sigma=maxDomain-minDomain;      %初始化当前尺度（标量）
optimalSolution=unifrnd(minDomain,maxDomain,DIM,groupNum); %unifrnd生成连续均匀分布的随机数，当前k个最优解的坐标（矩阵）
w=0;    %函数进化次数计数(标量)
for k=1:groupNum    
    funcV(k)=func(optimalSolution(:,k),DIM);
    w=w+1;
end      
[v_min,index_min]=min(funcV);  %找到k个点当中的最小值
while 1 % 尺度迭代开始
    stoNum=unifrnd(-sigma/2,sigma/2,DIM,groupNum); %产生尺度范围内的随机坐标,这些数都是以0为中心的
    optimalSolution=(repmat((optimalSolution(:,index_min))',groupNum,1))'+stoNum; %以position为中心，在原坐标上加上随机数
    %越界处理，并重新计算
    for k=1:groupNum  
        for d=1:DIM
            if optimalSolution(d,k)>=maxDomain
                %samplePos(d,k)=rand.*(maxDomain-minDomain);   
                optimalSolution(d,k)=maxDomain;
            end
            if optimalSolution(d,k)<=minDomain
                %samplePos(d,k)=rand.*(-maxDomain+minDomain);  
                optimalSolution(d,k)=minDomain;
            end               
        end
        funcV(k)=func(optimalSolution(:,k),DIM);
        w=w+1;
    end
    [v_min2,index_min2]=min(funcV);
    %fprintf('v_min=%e, v_min2=%e,sigma=%e,w=%d\n',v_min,v_min2,sigma,w);   
    %两次的最小值没有发生变化，(因为随机性的存在使得v_min2与v_min几乎没有相等的时候，根本不会去尺度收缩，所以选择二者的差值达到一定小，就认为稳定)认为达到了暂稳态
    if abs(v_min2-v_min)<gap
        if sigma<=sigmaMin
            break;
        end
        sigma=sigma/2.0; %尺度收缩
        fprintf('sigma=%e\n',sigma);
    elseif v_min2<v_min
       index_min=index_min2;
       fprintf('sigma=%e   v_min2=%e   ',sigma,v_min2);
       v_min=v_min2;   
       fprintf('v_min=%e\n',v_min);
    end
    
end 
tot_time=toc;
gbestV=v_min;
gfe=w;
fprintf('sigma=%e\n',sigma);
fprintf('gfe=%d,gbestV=%e,toc_time=%e\n',gfe,gbestV,tot_time);

% tot_time(rep) = toc; 
% gbestV(rep)=min(funcV);
% 
% [v_min,index_min]=min(funcV);
% position=optimalSolution(:,index_min);
% 
% gfe(rep)=w;
% fprintf('No. %d run, DIM=%d, Global minimum=%e. FE=%d, time=%d\n',rep,DIM,gbestV(rep),w,toc);
% %position
% if (gbestV(rep)>1e-5)
%     badspot=badspot+1.0;
% end    
% 
% fprintf('sigma=%e',sigma);
% fprintf('\n');
% fprintf('DIM=%d Repeat=%d sigmaMin=%0.1e MeanFE=%1.2e,Meantime=%1.2e\n',DIM,repeat,sigmaMin,mean(gfe),mean(tot_time));
% fprintf('MeanValue=%1.2e, BestValue=%1.2e, Std=%1.2e, \n',mean(gbestV),min(gbestV),std(gbestV));
% GDpercent=1-badspot/repeat;
% fprintf('GDpercent= %d BadSpot= %d w=%d \n\n',GDpercent,badspot,w); 
% xlabel('重复次数');
% ylabel('最小值');
% semilogy(gbestV); 
% output=gbestV;
% end


%  a=globalMin;
%   save a