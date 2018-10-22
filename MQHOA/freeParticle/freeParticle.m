% function freeParticle
% ���ܼ��ȶ��汾
clear all;         %������б���������ȫ�ֱ���
%clf;
hold on;
feature jit off;   %�ر�jit������
format long;       %�����������ʾ��ʽ��λ��
global DIM;        %����ȫ�ֱ���DIM������Ϊ��Դ�ļ����������������ã���Ч��ΧΪ�Ӷ��������λ�ÿ�ʼ����Դ�ļ�������

DIM=300;             %Ŀ�꺯����ά��
repeat=30;          %�ظ�����Ĵ���  
sigmaMin=0.000001;  %��С�߶�
groupNum=12;        %��˹��������Ŀ��kֵ
minDomain=-5;maxDomain=5; %Ŀ�꺯���Ķ��������½�  
gbestV=zeros(1,repeat);      %ÿ�μ�������ź���ֵ
gfe=zeros(1,repeat);         %ÿ�μ���ĺ�����������
tot_time = zeros(1,repeat);  %ÿ�μ����ʱ��
MFE=9000;  %����������

badspot=0;
for rep=1:repeat
    tic;           %��ʱ��ʼ      
    funcV=1./zeros(1,groupNum);        %�洢k��������ĺ���ֵ��������
    samplePos=zeros(DIM,groupNum);  %��ǰk�����������꣨����
    sigma=maxDomain-minDomain;      %��ʼ����ǰ�߶ȣ�������
    stdPre=zeros(DIM,1);            %��һ��ÿ��ά���ϲ�����ı�׼�������
    stdNow=zeros(DIM,1);            %��ǰÿ��ά���ϲ�����ı�׼�������
   
    %optimalSolution=unifrnd(minDomain,maxDomain,DIM,groupNum); %unifrnd�����������ȷֲ������������ǰk�����Ž�����꣨����
    %���Ƚ����еĲ������ɢ�İ����ڶ�����������
    optimalSolution=zeros(DIM,groupNum);
    gap=(maxDomain-minDomain)/(groupNum-2);
    optimalSolution(:,1)=minDomain;
    optimalSolution(:,groupNum)=maxDomain;
    for i=1:groupNum-2
        optimalSolution(:,i+1)=optimalSolution(:,i)+gap;
    end
    
    stdPre=std(optimalSolution,1,2); %std���������׼��.std(A,1,1)�����������׼��
    w=0;    %����������������(����)
    for k=1:groupNum    
        funcV(k)=func(optimalSolution(:,k),DIM);
        w=w+1;
%         if w>MFE
%             break;
%         end
    end      
%     fprintf('%f ',funcV);
%     fprintf('\n');
    while 1 % �߶ȵ�����ʼ 
        while 1  %г���ӵ�����ʼ
            while 1  %{�ܼ��ȶ�����������ʼ}
                for k=1:groupNum     %�˲���ȡ����ĵ�Ļ��ǲ���Ч�������һЩ�أ����ô˴������ܶȲ�����
                    %���Ȳ���---
                    % �����������ȷֲ��������        
                    %samplePos(:,k)=unifrnd(optimalSolution(:,k)-sigma/2,optimalSolution(:,k)+sigma/2,DIM,1);
                    %samplePos(:,k)=random('Uniform',optimalSolution(:,k)-sigma/2,optimalSolution(:,k)+sigma/2,DIM,1);
                    uniRand=random('Uniform',-sigma/2,sigma/2,DIM,1);
                    samplePos(:,k)=optimalSolution(:,k)+uniRand;
                    %Ҳ������a+(b-a).*rand(m,n)
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
                        funcV(k)=sampleValue;  %���µ�k�������ĺ���ֵ
                        optimalSolution(:,k)=samplePos(:,k); %������õ�k����
                    end
                end %for k ��������
                stdNow=std(optimalSolution,1,2); %�½��׼��
                c = stdPre-stdNow; %�������Ž��׼��Ĳ�ֵ
                stdPre = stdNow;   %��һ�α�׼�����Ϊ��ǰ��׼��
                if max(abs(c))<=sigma 
                    break;
                end
           end %�ܼ��ȶ���������
           meanPos=mean(optimalSolution,2);     %ȡ��ƽ������
           [v_max,index_max]=max(funcV);        %ȡ�����ֵ�����index_max
           optimalSolution(:,index_max)=meanPos;%��ƽ�������滻���ֵ��Ӧ����          
           funcV(index_max)=func(meanPos,DIM);
           w=w+1;
           if w>MFE                
               break;
           end
           stdPre=std(optimalSolution,1,2);     %�½��׼��
           if max(stdPre)<sigma    %�߶��½��о�
               break;
           end
        end %г���ӵ�������  
        if sigma<=sigmaMin %�����о� 
           break;
        end
        sigma=sigma/2.0;            %�߶��½�����  
    end % �߶ȵ�������
    tot_time(rep) = toc; 
    gbestV(rep)=min(funcV);
    
    [v_min,index_min]=min(funcV);
    position=optimalSolution(:,index_min);
    
    gfe(rep)=w;
%     fprintf('k�������������ĺ���ֵ�ֱ�Ϊ��\n');
%     fprintf('%e ',funcV);
%     fprintf('\n');
    fprintf('No. %d run, DIM=%d, Global minimum=%e. FE=%d, time=%d\n',rep,DIM,gbestV(rep),w,toc);
    %position
    if (gbestV(rep)>1e-5)
        badspot=badspot+1.0;
    end
    
 end % �ظ�����
fprintf('\n');
fprintf('DIM=%d Repeat=%d sigmaMin=%0.1e MeanFE=%1.2e,Meantime=%1.2e\n',DIM,repeat,sigmaMin,mean(gfe),mean(tot_time));
fprintf('MeanValue=%1.2e, BestValue=%1.2e, Std=%1.2e, \n',mean(gbestV),min(gbestV),std(gbestV));
GDpercent=1-badspot/repeat;
fprintf('GDpercent= %d BadSpot= %d w=%d \n\n',GDpercent,badspot,w); 
xlabel('�ظ�����');
ylabel('��Сֵ');
semilogy(gbestV); 
output=gbestV;
% end


%  a=globalMin;
%   save a