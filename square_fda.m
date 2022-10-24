%% FDAƽ�������(M^2-FDA)
clc;clear ;close;

%% ------FDA�״��������
j=sqrt(-1);
M=12; %������Ԫ��Ŀ
f0=2e9; %�ز�����Ƶ��
delta_f=3000; %Ƶ��ƫ��
Delta_f=delta_f*(1:M).^2;%����Ƶƫ����
c=3e8;        %����
lamda=c/f0;  %����
d=lamda/2;    %��Ԫ���
D=d*(0:M-1);  %���о�������
%Ru=c/delta_f;  %�����ģ������
theta=(-90:1:90)*pi/180; %�����Ƕ�����
R=linspace(0,3e5,1000); %������������
f=f0+(0:M-1)*delta_f; %��Ԫ��Ƶ�����������������ӣ�
R0 = 1e5;
theta0 = 30/180*pi;  %%����ָ��Ŀ��ĽǶȺ;���


%% -----ƽ��FDA��������ͼ
P1 = zeros(length(theta),length(R)); %��������ͼ

 for n = 1 : length(theta)
    for m = 1 : length(R)
         a1=exp(-j*2*pi/c*(Delta_f'*R(m)-D'*f0*sin(theta(n)))); %����ʸ��
         w1=exp(-j*2*pi/c*(Delta_f'*R0-D'*f0*sin(theta0))); 
        P1(n,m) =w1'*a1;
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
% %% -----ƽ��FDAʱ��Ƕ�ά��������ͼ     ��r0,theta0,t0=0.2ms������Ŀ��
% T=linspace(0,2e-3,1000);
% P2 = zeros(length(theta),length(T)); %��������ͼ
%  for n = 1 : length(theta)
%     for m = 1 : length(T)
%          a2=exp(-j*2*pi/c*(-Delta_f'*T(m)*c-D'*f0*sin(theta(n)))); %����ʸ��
%          w2=exp(-j*2*pi/c*(-Delta_f'*T(100)*c-D'*f0*sin(theta0))); 
%         P2(n,m) =w2'*a2;
%     end
%  end
% %% ��ͼ��ʱ��Ƕ�ά
% P2=P2';
% figure(2); 
% imagesc(theta*180/pi,T,abs(P2)/max(max(abs(P2)))); 
% xlabel('\theta^o'); ylabel('ʱ��/ms'); 
% axis tight; axis xy;
% title('ƽ��-FDAʱ��Ƕ�ά��������ͼ');
% colorbar;

