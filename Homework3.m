%%-------------------------------------------------------------------------
% 作者：   赵敏琨
% 日期：   2021年4月7日
% 说明：   求平面上100个点的凸多边形闭包及该多边形面积
% 自编代码，请勿PO上网，算法参考https://blog.csdn.net/viafcccy/article/details/87483567
%%-------------------------------------------------------------------------
%% 求平面上100个点的凸多边形闭包及该多边形面积
clc, clear, close all

n = 100;
% rng(3);
% 200*200方框正态分布
x = -100 + (100 - (-100))*randn(n,1);
y = -100 + (100 - (-100))*randn(n,1);
% 200*200方框随机分布
% x = -100 + (100 - (-100))*rand(n,1);
% y = -100 + (100 - (-100))*rand(n,1);

% ---------MATLAB自带函数实现
% dt=delaunayTriangulation(x,y);
% k = convexHull(dt);
% plot(x,y,'*');
% hold on;
% plot(x(k),y(k),'r');

% ---------自编代码（Graham扫描法）
% 点结构体赋值
for m = 1:n
    Point(m).x = x(m);
    Point(m).y = y(m);
end

[ysort,ysortindex] = sort(y);  %纵坐标排序,将纵坐标最小的作为起始点
xmove = x-x(ysortindex(1));
ymove = y-ysort(1); %将纵坐标最小的作为原点
subplot(2,2,1)
plot(xmove,ymove,'*');   %平移后的点图
title('平移后的点图')

for m = 1:n
    PointMove(m).x = xmove(m);
    PointMove(m).y = ymove(m);
    if ymove(m) == 0
        y0index = m;    %原点索引
    end
    PAarray(m) = Point_PolarAngle(PointMove(m));  %辐角序列
end
PAarray(y0index) = -1; %原点的辐角设为-1°，作为原点的标志（其他辐角均为非负）
[PAsort, PAindex] =sort(PAarray);   %辐角有小到大排序

subplot(2,2,2)
plot(x,y,'*'),hold on
plot(x(PAindex),y(PAindex),'r');     %按逆时针方向连线
title('按逆时针方向连线')

% 下面的循环求出了包围所有点的闭环折线
for k = 3:n-1
    SPAi = size(PAindex);
    if k+1 > SPAi(2)
        break;
    end
    A = Point_Angle(Point(PAindex(k-1)),Point(PAindex(k)),Point(PAindex(k+1)));
    if A == 1
        PAindex(k)=[];  %右转顺时针消除索引
        k = k - 1;  
    end 
end
PAindexCloseLine = [PAindex y0index];

subplot(2,2,3)
plot(x,y,'*'),hold on
plot(x(PAindexCloseLine),y(PAindexCloseLine),'r')     %闭环连线
title('闭环连线外轮廓折线')

% 下面的循环求出凸包
k = 2;
while 1
    SPAi = size(PAindex);
    if k+1 > SPAi(2)
        for m = 2:SJde(2)-1
            %轮询每轮仅能消除一个内凹，且对应判断值会变成0
            %此处再算判断值是为了校正
            J(m) = Point_Angle(Point(PAindex(m-1)),Point(PAindex(m)),Point(PAindex(m+1)));
            if J(m)==1 %有1（即有右边/顺时针）则继续消凹
                break;
            end
        end
        k = 2;
        if J == 0 %全0则跳出while loop
            break;
        end
    end
    A = Point_Angle(Point(PAindex(k-1)),Point(PAindex(k)),Point(PAindex(k+1)));
    J(k) = Point_Angle(Point(PAindex(k-1)),Point(PAindex(k)),Point(PAindex(k+1)));
    if A == 1
        PAindex(k)=[];  %右转顺时针消除索引 
        J(k)=[];
        k = k - 1;  
    end
    SJde = size(J);
    k = k + 1;
end

PAindex = [PAindex y0index];    %闭包

subplot(2,2,4)
plot(x,y,'*'),hold on
plot(x(PAindex),y(PAindex),'r')     %凸包连线
title('凸多边形闭包')

TheArea = polyarea(x(PAindex),y(PAindex));
% disp(['凸多边形的面积是',num2str(TheArea)])
figure
plot(x,y,'*'),hold on
plot(x(PAindex),y(PAindex),'r')    %凸包连线
title(['凸多边形闭包的面积为',num2str(TheArea)])


%% 内嵌函数
function D = Point_Distance(PointS, PointE)
D = sqrt((PointS.x-PointE.x)^2+(PointS.y-PointE.y)^2);
end

function PA = Point_PolarAngle(PointA)   %反正切求辐角
PA = rad2deg(atan(PointA.y/PointA.x));
if PA < 0
    PA = 180 + PA;
end
end

function Bool = Point_Angle(PointS, PointC, PointE)    %判断直线左右
Sidecur = [PointC.x PointC.y]-[PointS.x PointS.y];  %当前边的向量
Sideaft = [PointE.x PointE.y]-[PointC.x PointC.y];  %下一条边的向量
Crorslt = Sidecur(1)*Sideaft(2)-Sidecur(2)*Sideaft(1);
if Crorslt > 0  %根据叉积正负判断左/右（逆/顺时针）
    Bool = 0;   %左转
else
    Bool = 1;   %右转
end
end
