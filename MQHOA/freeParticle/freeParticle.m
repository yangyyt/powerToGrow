% function freeParticle
% 有能级稳定版本
clear all;         %清除所有变量，包括全局变量
%clf;
hold on;
feature jit off;   %关闭jit加速器
format long;       %控制命令窗口显示方式和位数
global DIM;        %定义全局变量DIM，可以为本源文件中其它函数所公用，有效范围为从定义变量的位置开始到本源文件结束。

DIM=300;             %目标函数的维度
repeat=30;          %重复计算的次数  
sigmaMin=0.000001;  %最小尺度
groupNum=12;        %高斯采样的数目，k值
minDomain=-5;maxDomain=5; %目标函数的定义域上下界  
gbestV=zeros(1,repeat);      %每次计算的最优函数值
gfe=zeros(1,repeat);         %每次计算的函数进化次数
tot_time = zeros(1,repeat);  %每次计算的时间
MFE=9000;  %最大迭代次数

badspot=0;
for rep=1:repeat
    tic;           %计时开始      
    funcV=1./zeros(1,groupNum);        %存储k个采样点的函数值（向量）
    samplePos=zeros(DIM,groupNum);  %当前k个采样点坐标（矩阵）
    sigma=maxDomain-minDomain;      %初始化当前尺度（标量）
    stdPre=zeros(DIM,1);            %上一次每个维度上采样点的标准差（向量）
    stdNow=zeros(DIM,1);            %当前每个维度上采样点的标准差（向量）
   
    %optimalSolution=unifrnd(minDomain,maxDomain,DIM,groupNum); %unifrnd生成连续均匀分布的随机数，当前k个最优解的坐标（矩阵）
    %首先将所有的采样点分散的安排在定义域区间上
    optimalSolution=zeros(DIM,groupNum);
    gap=(maxDomain-minDomain)/(groupNum-2);
    optimalSolution(:,1)=minDomain;
    optimalSolution(:,groupNum)=maxDomain;
    for i=1:groupNum-2
        optimalSolution(:,i+1)=optimalSolution(:,i)+gap;
    end
    
    stdPre=std(optimalSolution,1,2); %std按照行求标准差.std(A,1,1)：按照列求标准差
    w=0;    %函数进化次数计数(标量)
    for k=1:groupNum    
        funcV(k)=func(optimalSolution(:,k),DIM);
        w=w+1;
%         if w>MFE
%             break;
%         end
    end      
%     fprintf('%f ',funcV);
%     fprintf('\n');
    while 1 % 尺度迭代开始 
        while 1  %谐振子迭代开始
            while 1  %{能级稳定收敛迭代开始}
                for k=1:groupNum     %此步采取更多的点的话是不是效果会更好一些呢，觉得此处采样密度不够大
                    %均匀采样---
                    % 产生连续均匀分布的随机数        
                    %samplePos(:,k)=unifrnd(optimalSolution(:,k)-sigma/2,optimalSolution(:,k)+sigma/2,DIM,1);
                    %samplePos(:,k)=random('Uniform',optimalSolution(:,k)-sigma/2,optimalSolution(:,k)+sigma/2,DIM,1);
                    uniRand=random('Uniform',-sigma/2,sigma/2,DIM,1);
                    samplePos(:,k)=optimalSolution(:,k)+uniRand;
                    %也可以用a+(b-a).*rand(m,n)
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
                    sampleValue=func(samplePos(:,k),DIM);
                    w=w+1;
                    if w>MFE
                        break;
                    end
                    if sampleValue<funcV(k)     
                        funcV(k)=sampleValue;  %更新第k个样本的函数值
                        optimalSolution(:,k)=samplePos(:,k); %更新最好的k个点
                    end
                end %for k 采样结束
                stdNow=std(optimalSolution,1,2); %新解标准差
                c = stdPre-stdNow; %两次最优解标准差的差值
                stdPre = stdNow;   %上一次标准差更新为当前标准差
                if max(abs(c))<=sigma 
                    break;
                end
           end %能级稳定迭代结束
           meanPos=mean(optimalSolution,2);     %取得平均坐标
           [v_max,index_max]=max(funcV);        %取得最大值的序号index_max
           optimalSolution(:,index_max)=meanPos;%用平均坐标替换最大值对应坐标          
           funcV(index_max)=func(meanPos,DIM);
           w=w+1;
           if w>MFE                
               break;
           end
           stdPre=std(optimalSolution,1,2);     %新解标准差
           if max(stdPre)<sigma    %尺度下降判据
               break;
           end
        end %谐振子迭代结束  
        if sigma<=sigmaMin %精度判据 
           break;
        end
        sigma=sigma/2.0;            %尺度下降操作  
    end % 尺度迭代结束
    tot_time(rep) = toc; 
    gbestV(rep)=min(funcV);
    
    [v_min,index_min]=min(funcV);
    position=optimalSolution(:,index_min);
    
    gfe(rep)=w;
%     fprintf('k个采样点迭代后的函数值分别为：\n');
%     fprintf('%e ',funcV);
%     fprintf('\n');
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