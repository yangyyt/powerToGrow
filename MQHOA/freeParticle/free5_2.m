%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%���������ﵽ����������ֻҪ���������ܴﵽ2�ξ�����---------ʵ������ﵽ�������ε���������������Ƚ�С,����������³߶Ȼ���������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;         %������б���������ȫ�ֱ���
clf;
feature jit off;   %�ر�jit������
format long;       %�����������ʾ��ʽ��λ��
global DIM;        %����ȫ�ֱ���DIM������Ϊ��Դ�ļ����������������ã���Ч��ΧΪ�Ӷ��������λ�ÿ�ʼ����Դ�ļ�������
DIM=50;            %Ŀ�꺯����ά��
repeat=5;          %�ظ�����Ĵ���  
sigmaMin=0.000001;  %��С�߶�
gap=0.0001;          %��̬��ֵ
groupNum=10;        %��˹��������Ŀ��kֵ
minDomain=-5;maxDomain=5; %Ŀ�꺯���Ķ��������½�  
gbestV=zeros(1,repeat);      %ÿ�μ�������ź���ֵcl
gfe=zeros(1,repeat);         %ÿ�μ���ĺ�����������
sig=zeros(1,repeat);         %��¼ÿ�ε������������sigma
tot_time = zeros(1,repeat);  %ÿ�μ����ʱ��
MFE=10000000;  %����������
badspot=0;
for rep=1:repeat
    tic;           %��ʱ��ʼ  
    num=0;         %��¼�ﵽ���������Ĵ���
    flag1=1;       %��һ���Ƿ������ı�־
    funcV=1./zeros(1,groupNum);        %�洢k��������ĺ���ֵ��������
    samplePos=zeros(DIM,groupNum);  %��ǰk�����������꣨����
    sigma=maxDomain-minDomain;      %��ʼ����ǰ�߶ȣ�������
    optimalSolution=unifrnd(minDomain,maxDomain,DIM,groupNum); %unifrnd�����������ȷֲ������������ǰk�����Ž�����꣨����
    w=0;    %����������������(����)
    for k=1:groupNum     %����k��������ĺ���ֵ
        funcV(k)=func(optimalSolution(:,k),DIM);
        w=w+1;
    end      
    [v_min,index_min]=min(funcV);  %�ҵ�k���㵱�е���Сֵ
    while 1 % �߶ȵ�����ʼ
        for k=1:groupNum 
            %���ȷֲ�����
            stoNum=unifrnd(-sigma/2,sigma/2);
            optimalSolution(:,k)=optimalSolution(:,index_min)+stoNum;
            %Խ�紦���ͼ����µĲ�����ĺ���ֵ
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
        
        if sigma<=sigmaMin | w>MFE  %�߶��������㹻С���ߴﵽ�趨�����������������˳�
            break;
        end
        [v_min2,index_min2]=min(funcV);
        %���ε���Сֵû�з����仯��(��Ϊ����ԵĴ���ʹ��v_min2��v_min����û����ȵ�ʱ��Ҳ�ͻ�������߶�����������ѡ����ߵĲ�ֵ�ﵽһ��С������Ϊ�ȶ�)��Ϊ�ﵽ������̬
        if abs(v_min2-v_min)<=gap
            num=num+1;
            if flag1==1   %��ʾ�ϴδﵽ���������ж������Ҳ�ﵽ�������ж���num++���ﵽ2��sigma������������num���㣬Ϊ����ʵ����׼����
                if num==2
                    sigma=sigma/2.0; %�߶�����
                    num=0;
                end
            else         %��δﵽ�������ж����ϴ�û�ﵽ�����ж���flag1֮ǰ��0��������Ҫ���ı�־Ϊ1
                flag1=1; 
            end  
        else
            flag1=0;   %�˴�δ�ﵽ�������������ñ�־Ϊ0��ͬʱnum���㣬��ʾ�˴θ��ϴε�v_min���󣬴ﲻ������������
            num=0;
        end
        if v_min2<v_min  %����˴ε���СֵС��ǰһ�ε���Сֵ�����滻λ�ú�ֵ����Ϊ��һ�εĲ�������
           index_min=index_min2;
           v_min=v_min2;   
        end
    end 
    tot_time(rep) = toc; 
    gbestV(rep)=v_min;  %ÿһ�ε����ҵ�����Сֵ
    gfe(rep)=w;        %ÿһ�ε����ĺ�����������
    sig(rep)=sigma;    %ÿһ�ε����������������sigmaֵ
    fprintf('No. %d run, DIM=%d, Global minimum=%e. FE=%d, time=%d, sigma=%0.1e\n',rep,DIM,v_min,w,toc,sigma);
    if (gbestV(rep)>1e-5)
        badspot=badspot+1.0;
    end    
end
fprintf('�������ӣ����������������������ͳ߶�������Ҫ���������δﵽ����������Ϊ����');
fprintf('\n');
fprintf('DIM=%d Repeat=%d sigmaMin=%0.1e MeanFE=%1.2e,Meantime=%1.2e\n',DIM,repeat,sigmaMin,mean(gfe),mean(tot_time));
fprintf('MeanValue=%1.2e, BestValue=%1.2e, Std=%1.2e, \n',mean(gbestV),min(gbestV),std(gbestV));
GDpercent=1-badspot/repeat;
fprintf('GDpercent= %d BadSpot= %d w=%d \n\n',GDpercent,badspot,w); 
xlabel('�ظ�����');
ylabel('��Сֵ');
semilogy(gbestV); 
output=gbestV;

a=sig;
save a;
b=gbestV;
save b;