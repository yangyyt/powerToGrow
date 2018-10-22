% function freeParticle
%�߶Ƚ�����ȥ��û���߹�ȥ
clear all;         %������б���������ȫ�ֱ���
clf;
%hold on;
feature jit off;   %�ر�jit������
format long;       %�����������ʾ��ʽ��λ��
global DIM;        %����ȫ�ֱ���DIM������Ϊ��Դ�ļ����������������ã���Ч��ΧΪ�Ӷ��������λ�ÿ�ʼ����Դ�ļ�������
DIM=10;             %Ŀ�꺯����ά��
repeat=5;          %�ظ�����Ĵ���  
sigmaMin=0.0000001;  %��С�߶�
gap=0.00001;   %��̬��ֵ
groupNum=12;        %��˹��������Ŀ��kֵ
minDomain=-5;maxDomain=5; %Ŀ�꺯���Ķ��������½�  
gbestV=zeros(1,repeat);      %ÿ�μ�������ź���ֵ
gfe=zeros(1,repeat);         %ÿ�μ���ĺ�����������
sig=zeros(1,repeat);         %��¼ÿ�ε������������sigma
tot_time = zeros(1,repeat);  %ÿ�μ����ʱ��
MFE=10000000;  %����������
badspot=0;
for rep=1:repeat
    tic;           %��ʱ��ʼ  
    flag1=1;
    num=0;
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
    fprintf('v_min=%f   \n',v_min);
%     position=optimalSolution(:,index_min)
%     position 
%     fprintf('�����ʼ�����������£�\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for i=1:groupNum
%         for j=1:DIM
%             fprintf('%f ',optimalSolution(j,i));
%         end
%         fprintf('\n');
%     end
    while 1 % �߶ȵ�����ʼ
        %stoNum=unifrnd(-sigma/2,sigma/2,DIM,groupNum); %�����߶ȷ�Χ�ڵ��������,��Щ��������0Ϊ���ĵ�
        stoNum=unifrnd(-sigma/2,sigma/2);
        optimalSolution=(repmat((optimalSolution(:,index_min))',groupNum,1))'+stoNum; %��positionΪ���ģ���ԭ�����ϼ��������
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         fprintf('���������е��������£�\n');
%         for i=1:groupNum
%             for j=1:DIM
%                 fprintf('%f ',optimalSolution(j,i));
%             end
%             fprintf('\n');
%         end
        %Խ�紦�������¼���
        for k=1:groupNum  
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
        [v_min2,index_min2]=min(funcV);
%         if abs(v_min2-v_min)>gap
%             flag1=0;
%             num=0;
%         else
%             if flag1==1
%                 num=num+1;
%                 if num==2
%                     sigma=sigma/2.0; %�߶�����
%                 end
%             else
%                 flag1=1;
%                 num=num+1;
%             end
%             
%         end
        if abs(v_min2-v_min)<gap
            sigma=sigma/2.0;
        end
        if v_min2<v_min
           fprintf('-----------------------v_min2<v_min-------�����滻--------------------------------\n');
           fprintf('�˴ε�v_min�ǣ�%e  sigma=%e\n',v_min,sigma);
           index_min=index_min2;
           v_min=v_min2;   
           fprintf('�滻���v_min�ǣ�%e  ��ʱȡ����Сֵ�������ǣ�\n',v_min);
           for i=1:DIM
                fprintf('%f  ',optimalSolution(i,index_min));
           end
           fprintf('\n');
        end
%         position=optimalSolution(:,index_min2);
        
        fprintf('v_min2=%e    ',v_min2);
       
        %position
        %���ε���Сֵû�з����仯��(��Ϊ����ԵĴ���ʹ��v_min2��v_min����û����ȵ�ʱ�򣬸�������ȥ�߶�����������ѡ����ߵĲ�ֵ�ﵽһ��С������Ϊ�ȶ�)��Ϊ�ﵽ������̬
%         if abs(v_min2-v_min)<gap
%             if sigma<=sigmaMin
%                 break;
%             end
%             sigma=sigma/2.0; %�߶�����
%         end

        
       
%         if w>MFE
%             break;
%         end
    end 
    tot_time(rep) = toc; 
    gbestV(rep)=v_min;
%     [v_min,index_min]=min(funcV);
%     position=optimalSolution(:,index_min);
    gfe(rep)=w;
    sig(rep)=sigma;    %ÿһ�ε����������������sigmaֵ
    fprintf('No. %d run, DIM=%d, Global minimum=%e. FE=%d, time=%d, sigma=%0.1e \n',rep,DIM,v_min,w,toc,sigma);
    %position
    if (gbestV(rep)>1e-5)
        badspot=badspot+1.0;
    end    
end
fprintf('sigma=%e',sigma);
fprintf('\n');
fprintf('DIM=%d Repeat=%d sigmaMin=%0.1e MeanFE=%1.2e,Meantime=%1.2e\n',DIM,repeat,sigmaMin,mean(gfe),mean(tot_time));
fprintf('MeanValue=%1.2e, BestValue=%1.2e, Std=%1.2e, \n',mean(gbestV),min(gbestV),std(gbestV));
%fprintf('gbestV=%1.2e   ',gbestV);
GDpercent=1-badspot/repeat;
fprintf('GDpercent= %d BadSpot= %d w=%d \n\n',GDpercent,badspot,w); 
xlabel('�ظ�����');
ylabel('��Сֵ');
semilogy(gbestV); 
output=gbestV;
% end


%  a=globalMin;
%   save a