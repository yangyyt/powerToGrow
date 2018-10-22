% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Ackley function.
% % The global minima: f(x)=0  
% % Position (0)  mqhoaA300ά12�������㶼�ܴﵽ�ܸߵľ���
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = func(x,DIM)
% n = DIM;
% a = 20; b = 0.2; c = 2*pi;
% s1 = 0; s2 = 0;
% for i=1:n;
%    s1 = s1+(x(i))^2;
%    s2 = s2+cos(c*(x(i)));
% end
% y = -a*exp(-b*sqrt(1/n*s1))-exp(1/n*s2)+a+exp(1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Alpine function
% % The global minima: f(x)=0
% % Position (0)                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = func(x,DIM)
% n = DIM; 
% s=0;
%  for i = 1:DIM 
%         s = s + abs(x(i)*sin(x(i)) + 0.1*x(i));        
%  end
%  y=s;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Griewank function
% % The global minima: f(x)=0
% % Position (0)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = func(x,DIM)
% n = DIM;
% fr = 4000;
% s = 0;
% p = 1;
% for j = 1:n; s = s+(x(j))^2; end
% for j = 1:n; p = p*cos((x(j))/sqrt(j)); end
% y = s/fr-p+1;
% 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Sphere function 
% % The global minima: f(x)=0
% % Position (0)
% ��άʱmqhoa�Ĳ����Ƚϵ�ʱ��Ч������Խ�෴������  Ч��������
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = func(x,DIM)
% n = DIM;
% s = 0;
% for j = 1:n
%     s = s+(x(j))^2; 
% end
% y = s;
 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rastrigin function   
% The global minima: f(x)=0
% Position (0)   һά����պϣ�����2ά���Ͻ����ʼ���ã������������ڼ���άʱ������852���ܴﵽ���ȣ�
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = func(x,DIM)
% n = DIM; 
% s = 0;
% for j = 1:n
%     s = s+((x(j)-2)^2-10*cos(2*pi*(x(j)-2))); 
% end
% y = 10*n+s;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Dixon and Price function.
% % The global minima: f(x)=0   
% % Position (0)
% ֻ��1άЧ���ã�����2ά�ģ������ܶȼӴ�Ҳ�����ҵ��ܺõ�Ч��.���Ƕ���2ά��mqhoaA�ڲ����ܶ�8ʱ�ܴﵽ10-12,8000����ʱ10-1,freeҲһ�����ڲ�����ܴ�ʱ��Ч����������
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = func(x,DIM)
% n = DIM;
% s1 = 0;
% for j = 2:n;
%     s1 = s1+j*(2*x(j)^2-x(j-1))^2;    
% end
% y = s1+(x(1)-1)^2;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Zakharov function.
% % The global minima: f(x)=0
% % Position (0)            �����������ڼ���ά������Ҳ�㲻���ܺõľ���
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = func(x,DIM)
% n = DIM;
% s1 = 0;
% s2 = 0;
% for j = 1:n;
%     s1 = s1+(x(j)-1)^2;
%     s2 = s2+0.5*j*(x(j)-1);
% end
% y = s1+s2^2+s2^4;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Sum Squares Function 
% % The global minima: f(x)=0
% % Position (0)                                 Ч���ܺ�
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = func(x,DIM)
% n = DIM;
% s = 0;
% for j = 1:n  
%     s=s+j*x(j)^2; 
% end
% y = s;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Quadric function
% % The global minima: f(x)=0
% % Position (0)                   Ч��һ����
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = func(x,DIM)
% n = DIM; 
% s=0;
% s1 = 0;
%  for i = 1:DIM
%         for j=1:i
%              s1 = s1+x(j); 
%         end
%         s=s+s1^2;
%         s1=0;
%  end
%  y=s;  

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Levy function.
% % The global minima: f(x)=0
% % Position (1)
% ��άЧ�����ã�����5,6,7άʱmqhoaA�ɹ�100%��free0%????????���Ӳ����ܶ�mqhoa��200άʱ8000���������Ҳֻ�ܴﵽ10-5
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = func(x,DIM)
% n = DIM;
% for i = 1:n 
%     z(i) = 1+(x(i)-1-2)/4; 
% end
% s = sin(pi*z(1))^2;
% for i = 1:n-1
%     s = s+(z(i)-1)^2*(1+10*(sin(pi*z(i)+1))^2);
% end 
% y = s+(z(n)-1)^2*(1+(sin(2*pi*z(n)))^2);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Rosenbrock function 
% % The global minima: f(x)=0
% % Position (1)
% % DIM>1              
%groupNum     %kֵ10000�ǳ���׼8550Ҳ�ǳ���
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = func(x,DIM)
n = DIM;
sum = 0;
for j = 1:n-1;
    sum = sum+100*(x(j)^2-x(j+1))^2+(x(j)-1)^2;
end
y = sum;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Schwefel function
% % The global minima: f(x)=0
% % Position (420.9687)
% freeֻ��ʹ���ȴﵽ10��-4��mqhoaA:DIM=20;groupNum=8550;minDomain=-500;maxDomain=500;MFE=8500;
% ���Դﵽ100%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = func(x,DIM)
% n = DIM;
% sum = 0;
% for i=1:n
%     sum = sum+(x(i)*sin(sqrt(abs(x(i)))));
% end
% y = 418.9829*n-sum;

%��ά�����ֽ���
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Drop-Wave function
% % The global minima: f(x)=-1
% % Position (0)
% % [-5.12,5.12]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function  y = func(x,DIM)
% x1 = x(1);
% x2 = x(2);
% frac1 = 1 + cos(12*sqrt(x1^2+x2^2));
% frac2 = 0.5*(x1^2+x2^2) + 2;
% y = -frac1/frac2;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Beale function
% % The global minima: f(3,0.5)=0
% % Position (3,0.5)
% % DIM=2
% % [-4.5,4.5]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% free��΢��һЩ
% function  y = func(x,DIM)
% x1 = x(1);
% x2 = x(2);
% frac1 =(1.5-x1+x1*x2)^2 ;
% frac2 = (2.25-x1+x1*x2^2)^2;
% frac3=(2.625-x1+x1*x2^3)^2;
% y = frac1+frac2+frac3;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Three-hump camel  function
% % The global minima: f(0)=0
% % Position (0)
% % DIM=2
% % [-5,5]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function  y = func(x,DIM)
% x1 = x(1)-5;  %(5,5)�ڴ˵�ȡ����Сֵ0
% x2 = x(2)-5;
% y = 2*x1^2-1.05*x1^4+x1^6/6+x1*x2+x2^2;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % McCormick  function
% % The global minima: f(-0.54717,-1.54719)=-1.9133
% % Position (0)
% % DIM=2
% % [-4,4]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function  y = func(x,DIM)
% x1 = x(1);
% x2 = x(2);
% y = sin(x1+x2)+(x1-x2)^2-1.5*x1+2.5*x2+1;
