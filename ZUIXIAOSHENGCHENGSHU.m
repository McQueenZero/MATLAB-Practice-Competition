%% 最小生成树问题
clc,clear,close all
ZHANDIANImport = importdata('最小生成树坐标.txt');
Waterstation.num = ZHANDIANImport.textdata(:,1);
Waterstation.num(1) = [];
Waterstation.type = ZHANDIANImport.textdata(:,2);
Waterstation.type(1) = [];
Waterstation.x = ZHANDIANImport.data(:,1);
Waterstation.y = ZHANDIANImport.data(:,2);

% 点结构体赋值
for k = 1:13
    Point(k).num = k-1;
    Point(k).x = Waterstation.x(k);
    Point(k).y = Waterstation.y(k);
end

plot(Waterstation.x(1),Waterstation.y(1),'o'),hold on
plot(Waterstation.x(2:13),Waterstation.y(2:13),'*'),hold on
% for k = 1:13
%     text(Waterstation.x(k),Waterstation.y(k),num2str(k));
% end

DistanceMatrix = zeros(13,13);  %距离矩阵（代价矩阵）
for k = 2:13                %中心水站到一级水站距离
    DistanceMatrix(1,k) = Point_Distance(Point(1),Point(k));
end
for k = 2:13                %一级水站到中心水站设为无穷大
    DistanceMatrix(k,1) = Inf;
end
for k = 2:13
    for m = k+1:13
        DistanceMatrix(k,m) = Point_Distance(Point(k),Point(m));
        DistanceMatrix(m,k) = DistanceMatrix(k,m);
    end
end

[Tree,Dis] = MinSpanTree(DistanceMatrix,1);

Treesize = size(Tree);
k = 1;
for m = 1:Treesize(1)
    plot(Waterstation.x(Tree(k,:)),Waterstation.y(Tree(k,:)),'g','LineWidth',2),hold on
    k = k+1;
end

%% 内嵌函数
function D = Point_Distance(PointS, PointE)
D = sqrt((PointS.x-PointE.x)^2+(PointS.y-PointE.y)^2);
end

function [Tree,Dis]=MinSpanTree(adjMat,sPoint)
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

