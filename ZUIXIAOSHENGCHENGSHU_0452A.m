%% 最小生成树问题
clc,clear,close all
ZHANDIANImport = importdata('最小生成树坐标0452A.txt');
Waterstation.num = ZHANDIANImport.textdata(:,1);
Waterstation.type = ZHANDIANImport.textdata(:,2);
Waterstation.x = ZHANDIANImport.data(:,1);
Waterstation.y = ZHANDIANImport.data(:,2);

% 点结构体赋值
for k = 1:181
    Point(k).num = k-1;
    Point(k).x = Waterstation.x(k);
    Point(k).y = Waterstation.y(k);
end

plot(Waterstation.x(1),Waterstation.y(1),'ko'),hold on
plot(Waterstation.x(2:13),Waterstation.y(2:13),'r*'),hold on
plot(Waterstation.x(14:181),Waterstation.y(14:181),'b.'),hold on

% for k = 1:181
%     text(Waterstation.x(k),Waterstation.y(k),num2str(k));
% end

DistanceMatrix = zeros(181,181);  %距离矩阵
for k = 2:13
    DistanceMatrix(1,k) = Point_Distance(Point(1),Point(k));    %中心水站到一级水站距离
    DistanceMatrix(k,1) = Inf;  %一级水站到中心水站设为无穷大
    for m = 14:181
        DistanceMatrix(1,m) = Inf;  %中心水站到二级水站设为无穷大
        DistanceMatrix(k,m) = Point_Distance(Point(k),Point(m));    %计算一级到二级水站距离
        DistanceMatrix(m,k) = Inf;  %二级水站到一级水站设为无穷大
        DistanceMatrix(m,1) = Inf;  %二级水站到中心水站设为无穷大
    end
end
for k = 2:13
    for m = k+1:13
        DistanceMatrix(k,m) = Point_Distance(Point(k),Point(m));
        DistanceMatrix(m,k) = DistanceMatrix(k,m);
    end
end
for k = 14:181
    for m = k+1:181
        DistanceMatrix(k,m) = Point_Distance(Point(k),Point(m));
        DistanceMatrix(m,k) = DistanceMatrix(k,m);
    end
end

[Tree,Dis] = MinSpanTree(DistanceMatrix,1);
[Tree1,Dis1] = MinSpanTree(DistanceMatrix,1);     %单起点，用于一级供水站
[Tree2,Dis2] = MinSpanTree(DistanceMatrix,1:13);  %多起点，用于二级供水站
D = 0; %总里程

Treesize = size(Tree);
k = 1;
for m = 1:Treesize(1)
    if Tree1(k,:) <= 13
        plot(Waterstation.x(Tree1(k,:)),Waterstation.y(Tree1(k,:)),'r','LineWidth',1),hold on
        D = D + Dis1(k);
    end
    if (Tree2(k,1) <= 13 && Tree2(k,2) > 13) || (Tree2(k,2) <= 13 && Tree2(k,1) > 13)
        plot(Waterstation.x(Tree2(k,:)),Waterstation.y(Tree2(k,:)),'b','LineWidth',1),hold on
        D = D + Dis2(k);
    end
    if Tree2(k,:) > 13
        plot(Waterstation.x(Tree2(k,:)),Waterstation.y(Tree2(k,:)),'g','LineWidth',1),hold on
        D = D + Dis2(k);
    end
    k = k+1;
end

disp(['管道总里程为 ' num2str(D) ' 公里'])  

%% 内嵌函数
function D = Point_Distance(PointS, PointE)
D = sqrt((PointS.x-PointE.x)^2+(PointS.y-PointE.y)^2);
end

function [Tree,Dis]=MinSpanTree(adjMat,sPoint)      %最小生成树非常棒棒棒的算法！！！，既支持单起点，又支持多起点
Mat=adjMat;        %用作计算的邻接矩阵（为了不改变原矩阵）
Mat(Mat==0)=inf;   %将矩阵为0的值调整为inf
if nargin<=1
    MinNum=min(min(Mat));
    [S,~]=find(Mat==MinNum);
    S=S(1);
else
    S=sPoint;          %S指已经取到的节点集
                       %以sPoint为搜索起始点
end
N=size(adjMat,1);  %节点总数
Tree=zeros(N-1,2); %以sPoint为起点的最小生成树
Dis=zeros(N-1,1);  %最小生成树每一分支长度
for i=1:N-1
    Mat(:,S)=inf;                %将通向已在节点集的点的路径长度设置为inf（防止点被二次取到）
    tempMat=Mat(S,:);            %节点集S中的点到其他节点距离
    MinNum=min(min(tempMat));    %找到当前最短距离
    [m,n]=find(tempMat==MinNum); %找到当前最短距离对应是哪条路
    Tree(i,:)=[S(m(1)),n(1)];    %将最短路添加到树
    Dis(i)=adjMat(S(m(1)),n(1)); %将最短距离添加到距离
    S=[S,n(1)];                  %将新的节点添加到已取到节点
end
end

