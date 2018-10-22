% function mqhoaA1
%function output=mqhoaA1(k,func)
% 有能级稳定版本
clear all;         
clf;
%hold on;
feature jit off;   
format long;       %控制命令窗口显示方式和位数
global DIM;       

DIM=300;             %目标函数的维度
repeat=30;          %重复计算的次数  
sigmaMin=0.000001;  %最小尺度
groupNum=12;      %k值10000非常精准8550也非常好     200维852个采样点就能达到很好的精度
minDomain=-10;maxDomain=10; %目标函数的定义域上下界
gbestV=1./zeros(1,repeat);      %每次计算的最优函数值
gfe=zeros(1,repeat);         %每次计算的函数进化次数
tot_time = zeros(1,repeat);  %每次计算的时间
MFE=9000;

badspot=0;
for rep=1:repeat
    tic;           %计时开始      
    funcV=zeros(1,groupNum);        %存储k个采样点的函数值（向量）
    samplePos=zeros(DIM,groupNum);  %当前k个采样点坐标（矩阵）
    sigma=maxDomain-minDomain;      %当前尺度（标量）
    stdPre=zeros(DIM,1);            %上一次每个维度上采样点的标准差（向量）
    stdNow=zeros(DIM,1);            %当前每个维度上采样点的标准差（向量）
%     optimalSolution=unifrnd(minDomain,maxDomain,DIM,groupNum); %unifrnd生成连续均匀分布的随机数，当前k个最优解的坐标（矩阵）
    %首先将所有的采样点分散的安排在定义域区间上,实验证明这样的效果比上边的那一句初始化开始的点的效果好很多，精确度上去了
    optimalSolution=zeros(DIM,groupNum);
    gap=(maxDomain-minDomain)/(groupNum-2);
    optimalSolution(:,1)=minDomain;
    optimalSolution(:,groupNum)=maxDomain;
    for i=1:groupNum-2
        optimalSolution(:,i+1)=optimalSolution(:,i)+gap;
    end
    
    stdPre=std(optimalSolution,1,2) ; %std按照行求标准差.std(A,1,1)：按照列求标准差
    w=0;    %函数进化次数计数(标量)
    
    for k=1:groupNum 
        funcV(k)=func(optimalSolution(:,k),DIM);
        w=w+1;
    end       
    
    while 1 % 尺度迭代开始 
        while 1  %谐振子迭代开始
            while 1  %{能级稳定收敛迭代开始}
                for k=1:groupNum
                    %采用Box-Muller方法生成DIM维新的正态分布采样点
                    theat=2*pi*rand(DIM,1);%rand：在0~1之间，生成DIM行1列的随机数。
                    R=sqrt(-2.0*log(rand(DIM,1)));
                    gaosiRand=R.*cos(theat);  
                    samplePos(:,k)=optimalSolution(:,k)+sigma*gaosiRand;
                    for d=1:DIM
                        if samplePos(d,k)>maxDomain
                            %samplePos(d,k)=maxDomain-rand.*(maxDomain-minDomain);    
                            samplePos(d,k)=maxDomain;
                        end
                        if samplePos(d,k)<minDomain
                            %samplePos(d,k)=minDomain+rand.*(maxDomain-minDomain);  
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
                        optimalSolution(:,k)=samplePos(:,k); 
                    end
                end %for k
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
        end                        %谐振子迭代结束  
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

% end


%  a=globalMin;
%   save a