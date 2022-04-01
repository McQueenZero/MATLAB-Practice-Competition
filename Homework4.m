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

Point_u = Point(1); %初始化点集
Point_v = Point(2:13);  %点集v需要每歩更新
Point_t = Point;  %临时点集t，包括代价信息，需要实时更新
Side = [];   %边集
n = 1;

for k = 1:13
    Point_t(k).target = 0;
    Point_t(k).cost = DistanceMatrix(1,k);
    Point_t(k).du = 0;
    costTemp(k) = Point_t(k).cost;
end
[P,I] = sort(costTemp);
vsize = size(Point_v);
for m = 1:vsize(2)
    numTemp(m) = Point_v(m).num;
end
SideTemp = [I(1) I(2)];
Point_t(I(1)).du = Point_t(I(1)).du+1;
Point_t(I(2)).du = Point_t(I(2)).du+1;
Side = [Side SideTemp];
Point_u = [Point_u Point_v(numTemp==I(2)-1)]; %最小值索引I(1)是自己，I(2)是最近点
Point_v(numTemp==I(2)-1) = []; 
n = n+1;

while n<12+4
    for k = 1:13
        if DistanceMatrix(I(2),k) < Point_t(k).cost %代价小于当上一状态的代价则更新
            Point_t(k).cost = DistanceMatrix(I(2),k);
            Point_t(k).target = I(2)-1;
        end
        if k == I(1)  %上一处理点代价设为Inf
            Point_t(k).cost = Inf;
        end
        if Point_t(k).du >= 3
            Point_t(k).cost = Inf;
        end
        costTemp(k) = Point_t(k).cost;
    end
    [P,I] = sort(costTemp);
    vsize = size(Point_v);
    numTemp = [];
    for m = 1:vsize(2)
        numTemp(m) = Point_v(m).num;
    end
    SideTemp = [Point_t(I(2)).target+1 I(2)];
    Point_t(Point_t(I(2)).target+1).du = Point_t(I(1)).du+1;
    Point_t(I(2)).du = Point_t(I(2)).du+1;
    Side = [Side; SideTemp];
    Point_u = [Point_u Point_v(numTemp==I(2)-1)]; %最小值索引I(1)是自己，I(2)是最近点
    Point_v(numTemp==I(2)-1) = [];
    n = n+1;    
end
Side(12,:) = [];
Sidesize = size(Side);
k = 1;
for m = 1:Sidesize(1)
    plot(Waterstation.x(Side(k,:)),Waterstation.y(Side(k,:)),'g','LineWidth',2),hold on
    k = k+1;
end

%% 内嵌函数
function D = Point_Distance(PointS, PointE)
D = sqrt((PointS.x-PointE.x)^2+(PointS.y-PointE.y)^2);
end
