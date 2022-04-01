clc,clear,close all
%% 第1题
% 画星形线
t=0:0.01*pi:2*pi;
subplot(1,2,1);
plot((cos(t)).^3,(sin(t)).^3)
grid on, axis tight
legend('星形线','Location','southoutside')
xlabel('x'),ylabel('y');
title('星形线')
% 画马鞍面
subplot(1,2,2);
x=-1:0.1:1;
y=-1:0.1:1;
[X,Y]=meshgrid(x,y);
Z=X.^2-Y.^2;
surf(X,Y,Z)
grid on, axis tight
colormap jet
colorbar SouthOutside
xlabel('x'),ylabel('y'),zlabel('z');
title('马鞍面')

%% 第2题
x=-1:0.1:1;
y=-1:0.1:1;
[X,Y]=meshgrid(x,y);
Z=2*X.^2+Y.^2;
surf(X,Y,Z)
grid on, axis tight
colormap hsv
colorbar EastOutside
xlabel('x'),ylabel('y'),zlabel('z');
title('曲面z=2x^2+y^2')

%% 第3题 元胞自动机：森林火灾模型
n = 200;     %元胞矩阵大小
Plight = 0.000001;  %闪电几率
Pgrowth = 0.001;  %生长几率
UL = [n 1:n-1]; %上和左邻域
DR = [2:n 1];   %下和右邻域
veg = zeros(n,n);        %初始化
% The value of veg:
% empty == 0  
% burning == 1
% green == 2
imh = image(cat(3,veg,veg,veg));    %串联图层并获取图像句柄
% veg==?返回逻辑值矩阵，veg(,)==?返回逻辑值
for k = 1:100000
    sum = (veg(UL,:) == 1) + (veg(:,UL) == 1) + (veg(DR,:) == 1) + (veg(:,DR) == 1);
    %根据规则更新森林矩阵：树 = 树 - 着火的树 + 新生的树
    veg = 2 * (veg == 2) - ( (veg == 2) & (sum > 0 | (rand(n,n) < Plight)) ) + 2 * ( (veg == 0) & rand(n,n) < Pgrowth);
    imh.CData=cat(3, (veg == 1), (veg == 2), zeros(n)); %更新并串联三个图层
    %                   红色分量    绿色分量    蓝色分量
    drawnow

    %跑完循环或ctrl+C终止
end

