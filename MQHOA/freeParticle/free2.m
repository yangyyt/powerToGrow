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
groupNum=12;        %高斯采样的数目，k值
minDomain=-5;maxDomain=5; %目标函数的定义域上下界  
gbestV=zeros(1,repeat);      %每次计算的最优函数值
gfe=zeros(1,repeat);         %每次计算的函数进化次数
tot_time = zeros(1,repeat);  %每次计算的时间
MFE=800;  %最大迭代次数

badspot=0;
for rep=1:repeat
    a=0;
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
    %position=optimalSolution(:,index_min);  %k个点中的最小值的坐标
    while 1 % 尺度迭代开始
        %首先从k个点中找到最小值的点，然后在以这个点为中心进行采样，采k个点，继续找其中的最小值，若最小值的位置不变（或者最小值不变）则进行尺度收缩
        %此处设置多次迭代（比如10次），若这些迭代次数找到的值都不变，则认为达到了该算法所求的最小，进而进行尺度收缩
        stoNum=unifrnd(-sigma/2,sigma/2,DIM,groupNum); %产生尺度范围内的随机坐标
%         repmat(optimalSolution(:,index_min)',groupNum,1); 
        optimalSolutionTo2=(repmat(optimalSolution(:,index_min)',groupNum,1))'+stoNum; %以position为中心
        %越界处理
        for k=1:groupNum
            for d=1:DIM
                if samplePos(d,k)>=maxDomain
                     %samplePos(d,k)=rand.*(maxDomain-minDomain);   
                    samplePos(d,k)=maxDomain;
                end
                if samplePos(d,k)<=minDomain
                    %samplePos(d,k)=rand.*(-maxDomain+minDomain);  
                    samplePos(d,k)=minDomain;
                end               
            end
            funcV(k)=func(optimalSolution(:,k),DIM);
            w=w+1;
            if w>MFE
                break;
            end
        end
        [v_min2,index_min2]=min(funcV);
        if sigma<=sigmaMin
            break;
        end
        if a<10
            if v_min2==v_min
               sigma=sigma/2.0;
               a=a+1;
            elseif v_min2<v_min
                index_min=index_min2;
                optimalSolution(:,index_min)=optimalSolution(:,index_min2);
            end
        end
    end 
    tot_time(rep) = toc; 
    gbestV(rep)=min(funcV);
    
    [v_min,index_min]=min(funcV);
    position=optimalSolution(:,index_min);

    gfe(rep)=w;
    fprintf('No. %d run, DIM=%d, Global minimum=%e. FE=%d, time=%d\n',rep,DIM,gbestV(rep),w,toc);
    %position
    if (gbestV(rep)>1e-5)
        badspot=badspot+1.0;
    end    
 end % 重复计算
fprintf('\n');
fprintf('DIM=%d Repeat=%d sigmaMin=%0.1e MeanFE=%1.2e,Meantime=%1.2e\n',DIM,repeat,sigmaMin,mean(gfe),mean(tot_time));
fprintf('MeanValue=%1.2e, BestValue=%1.2e, Std=%1.2e, \n',mean(gbestV),min(gbestV),std(gbestV));
GDpercent=1-badspot/repeat;
fprintf('GDpercent= %d BadSpot= %d w=%d \n\n',GDpercent,badspot,w); 
xlabel('重复次数');
ylabel('最小值');
semilogy(gbestV); 
output=gbestV;
% end


%  a=globalMin;
%   save a