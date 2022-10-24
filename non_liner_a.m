function a = non_liner_a(g,R,theta)
%NL_F 非线性频偏导向矢量（或矩阵）
%   g为非线性频偏的非线性函数；R 为距离；theta 为角度
c=3e8;
f0=2e9; %载波中心频率
j=sqrt(-1);
lamda=c/f0;  %波长
d=lamda/2;    %阵元间距
[N,M]=size(g);
delta_f=3000; %频率偏移



if N==1
D=d*(0:M-1);  %阵列距离设置
Delta_f=g*delta_f;
a=exp(-j*2*pi/c*(Delta_f'*R-D'*f0*sin(theta)));
else
    D=ones(N,M);
   for i=1:N
    D(i,:)=(0:M-1)*d;
   end
   Delta_f=g*delta_f;
   a=exp(-j*2*pi/c*(Delta_f'*R-D'*f0*sin(theta)));;
end

end

