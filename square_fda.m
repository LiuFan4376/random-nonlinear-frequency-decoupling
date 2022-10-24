%% FDA平方解耦合(M^2-FDA)
clc;clear ;close;

%% ------FDA雷达参数设置
j=sqrt(-1);
M=12; %发射阵元数目
f0=2e9; %载波中心频率
delta_f=3000; %频率偏移
Delta_f=delta_f*(1:M).^2;%阵列频偏设置
c=3e8;        %光速
lamda=c/f0;  %波长
d=lamda/2;    %阵元间距
D=d*(0:M-1);  %阵列距离设置
%Ru=c/delta_f;  %最大无模糊距离
theta=(-90:1:90)*pi/180; %测量角度向量
R=linspace(0,3e5,1000); %测量距离向量
f=f0+(0:M-1)*delta_f; %阵元载频向量（均匀线性增加）
R0 = 1e5;
theta0 = 30/180*pi;  %%天线指向目标的角度和距离


%% -----平方FDA波束方向图
P1 = zeros(length(theta),length(R)); %波束方向图

 for n = 1 : length(theta)
    for m = 1 : length(R)
         a1=exp(-j*2*pi/c*(Delta_f'*R(m)-D'*f0*sin(theta(n)))); %导向矢量
         w1=exp(-j*2*pi/c*(Delta_f'*R0-D'*f0*sin(theta0))); 
        P1(n,m) =w1'*a1;
    end
 end
 
%% 画图
P1=P1';
figure(1); 
imagesc(theta*180/pi,R,abs(P1)/max(max(abs(P1)))); 
xlabel('\theta^o'); ylabel('R/m'); 
axis tight; axis xy;
title('');
colorbar;
% %% -----平方FDA时间角度维波束方向图     在r0,theta0,t0=0.2ms处存在目标
% T=linspace(0,2e-3,1000);
% P2 = zeros(length(theta),length(T)); %波束方向图
%  for n = 1 : length(theta)
%     for m = 1 : length(T)
%          a2=exp(-j*2*pi/c*(-Delta_f'*T(m)*c-D'*f0*sin(theta(n)))); %导向矢量
%          w2=exp(-j*2*pi/c*(-Delta_f'*T(100)*c-D'*f0*sin(theta0))); 
%         P2(n,m) =w2'*a2;
%     end
%  end
% %% 画图：时间角度维
% P2=P2';
% figure(2); 
% imagesc(theta*180/pi,T,abs(P2)/max(max(abs(P2)))); 
% xlabel('\theta^o'); ylabel('时间/ms'); 
% axis tight; axis xy;
% title('平方-FDA时间角度维波束方向图');
% colorbar;

