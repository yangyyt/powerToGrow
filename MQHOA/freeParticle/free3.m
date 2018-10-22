% function freeParticle
clear all;         %������б���������ȫ�ֱ���
clf;
%hold on;
feature jit off;   %�ر�jit������
format long;       %�����������ʾ��ʽ��λ��
global DIM;        %����ȫ�ֱ���DIM������Ϊ��Դ�ļ����������������ã���Ч��ΧΪ�Ӷ��������λ�ÿ�ʼ����Դ�ļ�������
DIM=1;             %Ŀ�꺯����ά��
repeat=30;          %�ظ�����Ĵ���  
sigmaMin=0.000001;  %��С�߶�
gap=0.00001;   %��̬��ֵ
groupNum=12;        %��˹��������Ŀ��kֵ
minDomain=-5;maxDomain=5; %Ŀ�꺯���Ķ��������½�  
gbestV=zeros(1,repeat);      %ÿ�μ�������ź���ֵ
gfe=zeros(1,repeat);         %ÿ�μ���ĺ�����������
tot_time = zeros(1,repeat);  %ÿ�μ����ʱ��
%MFE=10000;  %����������
badspot=0;
tic;           %��ʱ��ʼ  
funcV=1./zeros(1,groupNum);        %�洢k��������ĺ���ֵ��������
samplePos=zeros(DIM,groupNum);  %��ǰk�����������꣨����
sigma=maxDomain-minDomain;      %��ʼ����ǰ�߶ȣ�������
optimalSolution=unifrnd(minDomain,maxDomain,DIM,groupNum); %unifrnd�����������ȷֲ������������ǰk�����Ž�����꣨����
w=0;    %����������������(����)
for k=1:groupNum    
    funcV(k)=func(optimalSolution(:,k),DIM);
    w=w+1;
end      
[v_min,index_min]=min(funcV);  %�ҵ�k���㵱�е���Сֵ
while 1 % �߶ȵ�����ʼ
    stoNum=unifrnd(-sigma/2,sigma/2,DIM,groupNum); %�����߶ȷ�Χ�ڵ��������,��Щ��������0Ϊ���ĵ�
    optimalSolution=(repmat((optimalSolution(:,index_min))',groupNum,1))'+stoNum; %��positionΪ���ģ���ԭ�����ϼ��������
    %Խ�紦�������¼���
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
    %���ε���Сֵû�з����仯��(��Ϊ����ԵĴ���ʹ��v_min2��v_min����û����ȵ�ʱ�򣬸�������ȥ�߶�����������ѡ����ߵĲ�ֵ�ﵽһ��С������Ϊ�ȶ�)��Ϊ�ﵽ������̬
    if abs(v_min2-v_min)<gap
        if sigma<=sigmaMin
            break;
        end
        sigma=sigma/2.0; %�߶�����
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
% xlabel('�ظ�����');
% ylabel('��Сֵ');
% semilogy(gbestV); 
% output=gbestV;
% end


%  a=globalMin;
%   save a