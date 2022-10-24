%% FDA�������Ƶƫ�����(RF-log-FDA)
clc;clear ;close;

%% ------FDA�״��������
j=sqrt(-1);
M=12; %������Ԫ��Ŀ
peak_3db=M/sqrt(2);
f0=2e9; %�ز�����Ƶ��
delta_f=5e3; %Ƶ��ƫ��
   Delta_f=delta_f*log(ceil(M*rand(100,M)));%����Ƶƫ����sin((0:M-1)*pi/(M-1))
%    Delta_f=log(M)*delta_f*sin(ceil(M*rand(100,M))*pi/(M-1));%����Ƶƫ����
c=3e8;        %����
lamda=c/f0;  %����
d=lamda/2;    %��Ԫ���
D=d*(0:M-1);  %���о�������
Ru=c/delta_f;  %�����ģ������
theta=(-90:1:90)*pi/180; %�����Ƕ�����
R=linspace(0,3e5,1000); %������������
% f=f0+(0:M-1)*delta_f; %��Ԫ��Ƶ�����������������ӣ�
R0 = 1e5; theta0 = 30/180*pi;   %%����ָ��Ŀ��ĽǶȺ;���



%% -----RF_log_FDA��������ͼ
P1 = zeros(length(theta),length(R)); %��������ͼ

 for n = 1 : length(theta)
    for m = 1 : length(R)
         a1=exp(-j*2*pi/c*(Delta_f'*R(m)-D'*f0*sin(theta(n)))); %����ʸ��
         w1=exp(-j*2*pi/c*(Delta_f'*R0-D'*f0*sin(theta0))); 
          P1(n,m) =dot(a1,w1)*ones(100,1)/100;
    end
 end
 
%% ��ͼ
P1=abs(P1');
figure(1); 
imagesc(theta*180/pi,R,abs(P1)/max(max(abs(P1)))); 
xlabel('\theta^o'); ylabel('R/m'); 
axis tight; axis xy;
title('CF-log-FDA���߷��䷽��ͼ');
colorbar;
[C,h]= contour(theta*180/pi,R,P1,[peak_3db peak_3db]);

% D=P1(:,121);
% figure(2);
% plot(R,abs(D)/max(abs(D)));
% xlabel('R/m');  ylabel('��һ������'); 
% title('����ά����ͼ')


%% -----RF_log_FDAʱ��Ƕ�ά��������ͼ     ��r0,theta0,t0=0.2ms������Ŀ��
% T=linspace(0,2e-3,1000);
% P2 = zeros(length(theta),length(T)); %��������ͼ
%  for n = 1 : length(theta)
%     for m = 1 : length(T)
%          a2=exp(-j*2*pi/c*(-Delta_f'*T(m)*c-D'*f0*sin(theta(n)))); %����ʸ��
%          w2=exp(-j*2*pi/c*(-Delta_f'*T(100)*c-D'*f0*sin(theta0))); 
%          for i=1:100
%              P2(n,m)=w2(:,i)'*a2(:,i)+P2(n,m);
%          end
%         P2(n,m) =P2(n,m)/i;
%     end
%  end
% %% ��ͼ��ʱ��Ƕ�ά
% P2=P2';
% figure(2); 
% imagesc(theta*180/pi,T,abs(P2)/max(max(abs(P2)))); 
% xlabel('\theta^o'); ylabel('ʱ��/ms'); 
% axis tight; axis xy;
% title('RF-log-FDAʱ��Ƕ�ά��������ͼ');
% colorbar;



% %% -----RF-log-FDAʱ�����ά��������ͼ     
% P3 = zeros(length(R),length(T)); %��������ͼ
%  for n = 1 : length(R)
%     for m = 1 : length(T)
%          a3=exp(-j*2*pi/c*(-Delta_f'*T(m)*c+Delta_f'*R(n))); %����ʸ��
%          w3=exp(-j*2*pi/c*(-Delta_f'*T(100)*c+Delta_f'*R0)); 
%         P3(n,m) =dot(a3,w3)*ones(100,1)/100;
%     end
%  end
% %% ��ͼ��ʱ�����ά��������ͼ
% P3=P3';
% figure(3); 
% imagesc(R,T,abs(P3)/max(max(abs(P3)))); 
% xlabel('R/m'); ylabel('ʱ��/ms'); 
% axis tight; axis xy;
% title('RF-log-FDAʱ�����ά��������ͼ');
% colorbar;

% 



