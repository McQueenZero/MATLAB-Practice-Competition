%% ��С����������
clc,clear,close all
ZHANDIANImport = importdata('��С����������0452A.txt');
Waterstation.num = ZHANDIANImport.textdata(:,1);
Waterstation.type = ZHANDIANImport.textdata(:,2);
Waterstation.x = ZHANDIANImport.data(:,1);
Waterstation.y = ZHANDIANImport.data(:,2);

% ��ṹ�帳ֵ
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

DistanceMatrix = zeros(181,181);  %�������
for k = 2:13
    DistanceMatrix(1,k) = Point_Distance(Point(1),Point(k));    %����ˮվ��һ��ˮվ����
    DistanceMatrix(k,1) = Inf;  %һ��ˮվ������ˮվ��Ϊ�����
    for m = 14:181
        DistanceMatrix(1,m) = Inf;  %����ˮվ������ˮվ��Ϊ�����
        DistanceMatrix(k,m) = Point_Distance(Point(k),Point(m));    %����һ��������ˮվ����
        DistanceMatrix(m,k) = Inf;  %����ˮվ��һ��ˮվ��Ϊ�����
        DistanceMatrix(m,1) = Inf;  %����ˮվ������ˮվ��Ϊ�����
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
[Tree1,Dis1] = MinSpanTree(DistanceMatrix,1);     %����㣬����һ����ˮվ
[Tree2,Dis2] = MinSpanTree(DistanceMatrix,1:13);  %����㣬���ڶ�����ˮվ
D = 0; %�����

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

disp(['�ܵ������Ϊ ' num2str(D) ' ����'])  

%% ��Ƕ����
function D = Point_Distance(PointS, PointE)
D = sqrt((PointS.x-PointE.x)^2+(PointS.y-PointE.y)^2);
end

function [Tree,Dis]=MinSpanTree(adjMat,sPoint)      %��С�������ǳ����������㷨����������֧�ֵ���㣬��֧�ֶ����
Mat=adjMat;        %����������ڽӾ���Ϊ�˲��ı�ԭ����
Mat(Mat==0)=inf;   %������Ϊ0��ֵ����Ϊinf
if nargin<=1
    MinNum=min(min(Mat));
    [S,~]=find(Mat==MinNum);
    S=S(1);
else
    S=sPoint;          %Sָ�Ѿ�ȡ���Ľڵ㼯
                       %��sPointΪ������ʼ��
end
N=size(adjMat,1);  %�ڵ�����
Tree=zeros(N-1,2); %��sPointΪ������С������
Dis=zeros(N-1,1);  %��С������ÿһ��֧����
for i=1:N-1
    Mat(:,S)=inf;                %��ͨ�����ڽڵ㼯�ĵ��·����������Ϊinf����ֹ�㱻����ȡ����
    tempMat=Mat(S,:);            %�ڵ㼯S�еĵ㵽�����ڵ����
    MinNum=min(min(tempMat));    %�ҵ���ǰ��̾���
    [m,n]=find(tempMat==MinNum); %�ҵ���ǰ��̾����Ӧ������·
    Tree(i,:)=[S(m(1)),n(1)];    %�����·��ӵ���
    Dis(i)=adjMat(S(m(1)),n(1)); %����̾�����ӵ�����
    S=[S,n(1)];                  %���µĽڵ���ӵ���ȡ���ڵ�
end
end

