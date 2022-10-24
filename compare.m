%% FDA对数log解耦合(log-FDA)
clc;clear ;close;

%% ------FDA雷达参数设置
j=sqrt(-1);
M=12; %发射阵元数目
peak_3db=M/sqrt(2);
f0=2e9; %载波中心频率
delta_f=5000; %频率偏移
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


%% -----log_FDA解耦波束方向图
g1=log((1:M));
P1 = zeros(length(theta),length(R)); %波束方向图
 for n = 1 : length(theta)
    for m = 1 : length(R)
         a1= non_liner_a(g1,R(m),theta(n)); %导向矢量
         w1= non_liner_a(g1,R0,theta0); 
        P1(n,m) =w1'*a1;
    end
 end
 
% %% 画图
 P1=abs(P1');
% figure(1)
% [C1,h]= contour(theta*180/pi,R,P1,[peak_3db peak_3db],'b','LineWidth',1);
% xlabel('\theta/^o'); ylabel('R/m');
% axis([-90,90,0,3e5]);
figure(1); 
imagesc(theta*180/pi,R,abs(P1)/max(max(abs(P1)))); 
xlabel('\theta^o'); ylabel('R/m'); 
axis tight; axis xy;
title('');
colorbar;



% %% -----sin_FDA解耦波束方向图
% g2=sin((0:M-1)*pi/(M-1));
% P2 = zeros(length(theta),length(R)); %波束方向图
% 
%  for n = 1 : length(theta)
%     for m = 1 : length(R)
%          a1=non_liner_a(g2,R(m),theta(n)); %导向矢量
%          w1=non_liner_a(g2,R0, theta0); 
%         P2(n,m) =w1'*a1;
%     end
%  end
%  
% %% 画图
% P2=P2';
% figure(2); 
% imagesc(theta*180/pi,R,abs(P2)/max(max(abs(P2)))); 
% xlabel('\theta^o'); ylabel('R/m'); 
% axis tight; axis xy;
% title('');
% colorbar;
% 
% %% -----平方FDA波束方向图
% P3 = zeros(length(theta),length(R)); %波束方向图
% g3=(1:M).^2;
%  for n = 1 : length(theta)
%     for m = 1 : length(R)
%         a1=non_liner_a(g3,R(m),theta(n)); %导向矢量
%          w1=non_liner_a(g3,R0, theta0); 
%         P3(n,m) =w1'*a1;
%     end
%  end
%  %% 画图
% P3=P3';
% figure(3); 
% imagesc(theta*180/pi,R,abs(P3)/max(max(abs(P3)))); 
% xlabel('\theta^o'); ylabel('R/m'); 
% axis tight; axis xy;
% title('');
% colorbar;
% -----RF_FDA波束方向图
     Delta_f=delta_f*ceil((M*rand(100,M)));%阵列频偏设置
     P4= zeros(length(theta),length(R)); %波束方向图

 for n = 1 : length(theta)
    for m = 1 : length(R)
         a1=exp(-j*2*pi/c*(Delta_f'*R(m)-D'*f0*sin(theta(n)))); %导向矢量
         w1=exp(-j*2*pi/c*(Delta_f'*R0-D'*f0*sin(theta0))); 
          P4(n,m) =dot(a1,w1)*ones(100,1)/100;
    end
 end
 
% %% 画图
% P3=abs(P3');
% % figure(3)
% % [C3,h]= contour(theta*180/pi,R,P3,[peak_3db peak_3db],'c','LineWidth',1);
% % xlabel('\theta/^o'); ylabel('R/m');
% % axis([-90,90,0,3e5]);
 %% 画图
P4=P4';
figure(4); 
imagesc(theta*180/pi,R,abs(P4)/max(max(abs(P4)))); 
xlabel('\theta^o'); ylabel('R/m'); 
axis tight; axis xy;
title('');
colorbar;

 %% -----RF_log_FDA波束方向图
     Delta_f=delta_f*log(ceil(M*rand(100,M)));%阵列频偏设置sin((0:M-1)*pi/(M-1))
     P5 = zeros(length(theta),length(R)); %波束方向图

 for n = 1 : length(theta)
    for m = 1 : length(R)
         a1=exp(-j*2*pi/c*(Delta_f'*R(m)-D'*f0*sin(theta(n)))); %导向矢量
         w1=exp(-j*2*pi/c*(Delta_f'*R0-D'*f0*sin(theta0))); 
          P5(n,m) =dot(a1,w1)*ones(100,1)/100;
    end
 end
 
%% 画图
P5=abs(P5');
% figure(4)
% [C2,h]= contour(theta*180/pi,R,P4,[peak_3db peak_3db],'r','LineWidth',1);
% xlabel('\theta/^o'); ylabel('R/m');
% axis([-90,90,0,3e5]);
figure(5); 
imagesc(theta*180/pi,R,abs(P5)/max(max(abs(P5)))); 
xlabel('\theta^o'); ylabel('R/m'); 
axis tight; axis xy;
% title('CF-log-FDA天线发射方向图');
colorbar;

% D=P4(:,121);
% figure(2);
% plot(R,abs(D)/max(abs(D)));
% xlabel('R/m');  ylabel('归一化幅度'); 
% title('距离维剖面图')

figure(6);
%  plot(R,abs(P1(:,121))/max(abs(P1(:,121))),'.-',...
%      R,abs(P2(:,121))/max(abs(P2(:,121))),'x-',...
%      R,abs(P3(:,121))/max(abs(P3(:,121))), 's-','LineWidth',1);
 plot(R,abs(P1(:,121))/max(abs(P1(:,121))),'.-',...
     R,abs(P4(:,121))/max(abs(P4(:,121))),'x-',...
     R,abs(P5(:,121))/max(abs(P5(:,121))), 's-','LineWidth',1)
xlabel('距离\m'); ylabel('归一化幅度'); 
% legend('对数频偏','正弦频偏','平方频偏','随机频偏','随机非线性频偏');
legend('对数频偏','随机频偏','随机非线性频偏');
axis([0,3e5,0,1]);
title('');
% axis([0.5e5,1.5e5,0,1]);
%  legend('对数频偏','正弦频偏','平方频偏');
% figure(6)
% plot(C1(1,2:end),C1(2,2:end),C2(1,2:end),C2(2,2:end),'r','LineWidth',1);
% xlabel('角度/^o'); ylabel('R/m');
% axis([-90,90,0,3e5]);
% legend('对数频偏','随机频偏','随机非线性频偏');




