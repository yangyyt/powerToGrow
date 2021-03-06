%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%只要达到一次收敛条件就尺度缩减
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;         %清除所有变量，包括全局变量
clf;
%hold on;
feature jit off;   %关闭jit加速器
format long;       %控制命令窗口显示方式和位数
global DIM;        %定义全局变量DIM，可以为本源文件中其它函数所公用，有效范围为从定义变量的位置开始到本源文件结束。
DIM=50;            %目标函数的维度
repeat=5;          %重复计算的次数  
sigmaMin=0.000001;  %最小尺度
gap=0.0001;          %稳态差值
groupNum=10;        %高斯采样的数目，k值
minDomain=-5;maxDomain=5; %目标函数的定义域上下界  
gbestV=zeros(1,repeat);      %每次计算的最优函数值
gfe=zeros(1,repeat);         %每次计算的函数进化次数
sig=zeros(1,repeat);         %记录每次迭代最后收敛的sigma
tot_time = zeros(1,repeat);  %每次计算的时间
MFE=10000000;  %最大迭代次数
badspot=0;
for rep=1:repeat
    tic;           %计时开始  
    num=0;         %记录达到收敛条件的次数
    flag1=1;       %上一次是否收敛的标志
    funcV=1./zeros(1,groupNum);        %存储k个采样点的函数值（向量）
    samplePos=zeros(DIM,groupNum);  %当前k个采样点坐标（矩阵）
    sigma=maxDomain-minDomain;      %初始化当前尺度（标量）
    optimalSolution=unifrnd(minDomain,maxDomain,DIM,groupNum); %unifrnd生成连续均匀分布的随机数，当前k个最优解的坐标（矩阵）
    w=0;    %函数进化次数计数(标量)
    for k=1:groupNum     %
        funcV(k)=func(optimalSolution(:,k),DIM);
        w=w+1;
    end      
    [v_min,index_min]=min(funcV);  %找到k个点当中的最小值
    while 1 % 尺度迭代开始
        for k=1:groupNum 
            %均匀分布采样
            stoNum=unifrnd(-sigma/2,sigma/2);
            optimalSolution(:,k)=optimalSolution(:,index_min)+stoNum;
            %越界处理，并计算新的采样点的函数值
            for d=1:DIM
                if optimalSolution(d,k)>=maxDomain 
                    optimalSolution(d,k)=maxDomain;
                end
                if optimalSolution(d,k)<=minDomain
                    optimalSolution(d,k)=minDomain;
                end               
            end
            funcV(k)=func(optimalSolution(:,k),DIM);
            w=w+1;
        end
        if sigma<=sigmaMin | w>MFE %尺度收缩到足够小或者迭代次数超过设定的最大迭代次数
            break;
        end
        %两次的最小值没有发生变化，(因为随机性的存在使得v_min2与v_min几乎没有相等的时候，根本不会尺度收缩，所以选择二者的差值达到一定小，就认为稳定)认为达到了暂稳态
        [v_min2,index_min2]=min(funcV);
        if abs(v_min2-v_min)<=gap
            sigma=sigma/2.0;
        end
        if v_min2<v_min  %如果此次的最小值小于前一次的最小值，则替换
           index_min=index_min2;
           v_min=v_min2;   
        end
    end 
    tot_time(rep) = toc; 
    gbestV(rep)=v_min;  %每一次迭代找到的最小值
    gfe(rep)=w;        %每一次迭代的函数迭代次数
    sig(rep)=sigma;    %每一次迭代的最后收敛到的sigma值
    fprintf('No. %d run, DIM=%d, Global minimum=%e. FE=%d, time=%d, sigma=%0.1e\n',rep,DIM,v_min,w,toc,sigma);
    if (gbestV(rep)>1e-5)
        badspot=badspot+1.0;
    end    
end
fprintf('sigma=%e',sigma);
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
a=sig;
save a;
b=gbestV;
save b;


%  a=globalMin;
%   save a