%% FDA���Ƶƫ�����(RF-FDA)
clc;clear ;close;

%% ------FDA�״��������
j=sqrt(-1);
M=12; %������Ԫ��Ŀ
f0=2e9; %�ز�����Ƶ��
delta_f=3000; %Ƶ��ƫ��
g=ceil((M*rand(500,M)));%����Ƶƫ����
c=3e8;        %����
lamda=c/f0;  %����
d=lamda/2;    %��Ԫ���
D=d*(0:M-1);  %���о�������
Ru=c/delta_f;  %�����ģ������
theta=(0:1:90)*pi/180; %�����Ƕ�����
R=linspace(0,3e5,1000); %������������
f=f0+(0:M-1)*delta_f; %��Ԫ��Ƶ�����������������ӣ�
R0 = 1e5;
theta0 = 30/180*pi;  %%����ָ��Ŀ��ĽǶȺ;���


%% -----CF_FDA���������ͼ

P1 = zeros(length(theta),length(R)); %��������ͼ

for n = 1 : length(theta)
    for m = 1 : length(R)
         a1=non_liner_a(g,R(m),theta(n)); %����ʸ��
         w1=non_liner_a(g,R0,theta0); 
          P1(n,m) =dot(a1,w1)*ones(500,1);
    end
 end
 
%% ��ͼ
P1=P1';
figure(1); 
imagesc(theta*180/pi,R,abs(P1)/max(max(abs(P1)))); 
xlabel('\theta^o'); ylabel('R/m'); 
axis tight; axis xy;
title('');
colorbar;


% %% -----RF-FDAʱ��Ƕ�ά��������ͼ     ��r0,theta0,t0=0.2ms������Ŀ��
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
% title('RF-FDAʱ��Ƕ�ά��������ͼ');
% colorbar;




% %% -----RF-FDAʱ�����ά��������ͼ     ��r0,theta0,t0=0.2ms������Ŀ��
% T=linspace(0,2e-3,1000);
% P3 = zeros(length(R),length(T)); %��������ͼ
%  for n = 1 : length(R)
%     for m = 1 : length(T)
%          a3=exp(-j*2*pi/c*(-Delta_f'*T(m)*c+Delta_f'*R(n))); %����ʸ��
%          w3=exp(-j*2*pi/c*(-Delta_f'*T(100)*c+Delta_f'*R0)); 
%         P3(n,m) =dot(a3,w3)*ones(100,1)/100;
%     end
%  end
% %% ��ͼ��ʱ��Ƕ�ά
% P3=P3';
% figure(3); 
% imagesc(R,T,abs(P3)/max(max(abs(P3)))); 
% xlabel('R/m'); ylabel('ʱ��/ms'); 
% axis tight; axis xy;
% title('RF-FDAʱ�����ά��������ͼ');
% colorbar;











