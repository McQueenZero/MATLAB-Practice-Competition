%% 一、Matlab数值运算
clc, clear, close all
A = magic(50);
I = find(A <= 100);
Acolsum = sum(A);   %列总和
Aallsum = sum(A(:));    %所有和
Alt100mean = mean(A(I));    %满足问题(2)平均值
A(I) = Alt100mean;

%% 二、Matlab数据可视化操作
clc, clear, close all
x = -pi/2:0.01*pi:pi/2;
nd = -pi/2 + pi*rand(1,10);
xd = -pi/2:pi/9:pi/2;

subplot(2,2,2)
plot(x,sqrt(cos(x)),'k'),hold on
plot(xd,nd,'r.')
title('y=\surd{cosx}')
axis tight

subplot(2,2,3)
t = 0:0.01*pi:2*pi;
xx = 2*cos(t);
yy = 4*sin(t);
plot(xx,yy)
axis([-5 5 -5 5])
title('x^2/2^2+y^2/4^2=1')

%% 三、Matlab三维绘图
clc, clear, close all
x = -9:0.5:9;
y = -9:0.5:9;
[X,Y] = meshgrid(x,y);
Z = (2*sin(X).*sin(Y))./(X.*Y);
% ZI = find(isnan(Z));
% Z(ZI) = 0;
mesh(X,Y,Z)

%% 四、Matlab符号计算
clc, clear, close all
syms x y
f1 = 1/x^2-3*x+3;
fcn1 = matlabFunction(f1);
df = diff(fcn1,x,3);
disp('函数的3阶导数是：')
disp(df)

f2 = (sin(x)).^2.*(cos(x)).^2;
fcn2 = matlabFunction(f2);
itgl = int(fcn2,x,-pi,pi);
itglvalue = eval(itgl);
disp('定积分的值是：')
disp(itglvalue)

f3 = (sin(x))^2+(sin(y))^2;
fcn3 = matlabFunction(f3);
mitgl = int(int(fcn3,y,1,x^2),x,1,2);
mitglvalue = eval(mitgl);
disp('二重积分的值是：')
disp(mitglvalue)

