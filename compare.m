%% FDA����log�����(log-FDA)
clc;clear ;close;

%% ------FDA�״��������
j=sqrt(-1);
M=12; %������Ԫ��Ŀ
peak_3db=M/sqrt(2);
f0=2e9; %�ز�����Ƶ��
delta_f=5000; %Ƶ��ƫ��
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


%% -----log_FDA���������ͼ
g1=log((1:M));
P1 = zeros(length(theta),length(R)); %��������ͼ
 for n = 1 : length(theta)
    for m = 1 : length(R)
         a1= non_liner_a(g1,R(m),theta(n)); %����ʸ��
         w1= non_liner_a(g1,R0,theta0); 
        P1(n,m) =w1'*a1;
    end
 end
 
% %% ��ͼ
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



% %% -----sin_FDA���������ͼ
% g2=sin((0:M-1)*pi/(M-1));
% P2 = zeros(length(theta),length(R)); %��������ͼ
% 
%  for n = 1 : length(theta)
%     for m = 1 : length(R)
%          a1=non_liner_a(g2,R(m),theta(n)); %����ʸ��
%          w1=non_liner_a(g2,R0, theta0); 
%         P2(n,m) =w1'*a1;
%     end
%  end
%  
% %% ��ͼ
% P2=P2';
% figure(2); 
% imagesc(theta*180/pi,R,abs(P2)/max(max(abs(P2)))); 
% xlabel('\theta^o'); ylabel('R/m'); 
% axis tight; axis xy;
% title('');
% colorbar;
% 
% %% -----ƽ��FDA��������ͼ
% P3 = zeros(length(theta),length(R)); %��������ͼ
% g3=(1:M).^2;
%  for n = 1 : length(theta)
%     for m = 1 : length(R)
%         a1=non_liner_a(g3,R(m),theta(n)); %����ʸ��
%          w1=non_liner_a(g3,R0, theta0); 
%         P3(n,m) =w1'*a1;
%     end
%  end
%  %% ��ͼ
% P3=P3';
% figure(3); 
% imagesc(theta*180/pi,R,abs(P3)/max(max(abs(P3)))); 
% xlabel('\theta^o'); ylabel('R/m'); 
% axis tight; axis xy;
% title('');
% colorbar;
% -----RF_FDA��������ͼ
     Delta_f=delta_f*ceil((M*rand(100,M)));%����Ƶƫ����
     P4= zeros(length(theta),length(R)); %��������ͼ

 for n = 1 : length(theta)
    for m = 1 : length(R)
         a1=exp(-j*2*pi/c*(Delta_f'*R(m)-D'*f0*sin(theta(n)))); %����ʸ��
         w1=exp(-j*2*pi/c*(Delta_f'*R0-D'*f0*sin(theta0))); 
          P4(n,m) =dot(a1,w1)*ones(100,1)/100;
    end
 end
 
% %% ��ͼ
% P3=abs(P3');
% % figure(3)
% % [C3,h]= contour(theta*180/pi,R,P3,[peak_3db peak_3db],'c','LineWidth',1);
% % xlabel('\theta/^o'); ylabel('R/m');
% % axis([-90,90,0,3e5]);
 %% ��ͼ
P4=P4';
figure(4); 
imagesc(theta*180/pi,R,abs(P4)/max(max(abs(P4)))); 
xlabel('\theta^o'); ylabel('R/m'); 
axis tight; axis xy;
title('');
colorbar;

 %% -----RF_log_FDA��������ͼ
     Delta_f=delta_f*log(ceil(M*rand(100,M)));%����Ƶƫ����sin((0:M-1)*pi/(M-1))
     P5 = zeros(length(theta),length(R)); %��������ͼ

 for n = 1 : length(theta)
    for m = 1 : length(R)
         a1=exp(-j*2*pi/c*(Delta_f'*R(m)-D'*f0*sin(theta(n)))); %����ʸ��
         w1=exp(-j*2*pi/c*(Delta_f'*R0-D'*f0*sin(theta0))); 
          P5(n,m) =dot(a1,w1)*ones(100,1)/100;
    end
 end
 
%% ��ͼ
P5=abs(P5');
% figure(4)
% [C2,h]= contour(theta*180/pi,R,P4,[peak_3db peak_3db],'r','LineWidth',1);
% xlabel('\theta/^o'); ylabel('R/m');
% axis([-90,90,0,3e5]);
figure(5); 
imagesc(theta*180/pi,R,abs(P5)/max(max(abs(P5)))); 
xlabel('\theta^o'); ylabel('R/m'); 
axis tight; axis xy;
% title('CF-log-FDA���߷��䷽��ͼ');
colorbar;

% D=P4(:,121);
% figure(2);
% plot(R,abs(D)/max(abs(D)));
% xlabel('R/m');  ylabel('��һ������'); 
% title('����ά����ͼ')

figure(6);
%  plot(R,abs(P1(:,121))/max(abs(P1(:,121))),'.-',...
%      R,abs(P2(:,121))/max(abs(P2(:,121))),'x-',...
%      R,abs(P3(:,121))/max(abs(P3(:,121))), 's-','LineWidth',1);
 plot(R,abs(P1(:,121))/max(abs(P1(:,121))),'.-',...
     R,abs(P4(:,121))/max(abs(P4(:,121))),'x-',...
     R,abs(P5(:,121))/max(abs(P5(:,121))), 's-','LineWidth',1)
xlabel('����\m'); ylabel('��һ������'); 
% legend('����Ƶƫ','����Ƶƫ','ƽ��Ƶƫ','���Ƶƫ','���������Ƶƫ');
legend('����Ƶƫ','���Ƶƫ','���������Ƶƫ');
axis([0,3e5,0,1]);
title('');
% axis([0.5e5,1.5e5,0,1]);
%  legend('����Ƶƫ','����Ƶƫ','ƽ��Ƶƫ');
% figure(6)
% plot(C1(1,2:end),C1(2,2:end),C2(1,2:end),C2(2,2:end),'r','LineWidth',1);
% xlabel('�Ƕ�/^o'); ylabel('R/m');
% axis([-90,90,0,3e5]);
% legend('����Ƶƫ','���Ƶƫ','���������Ƶƫ');




