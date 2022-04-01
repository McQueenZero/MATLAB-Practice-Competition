%% ��С����������
clc,clear,close all
ZHANDIANImport = importdata('��С����������.txt');
Waterstation.num = ZHANDIANImport.textdata(:,1);
Waterstation.num(1) = [];
Waterstation.type = ZHANDIANImport.textdata(:,2);
Waterstation.type(1) = [];
Waterstation.x = ZHANDIANImport.data(:,1);
Waterstation.y = ZHANDIANImport.data(:,2);

% ��ṹ�帳ֵ
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

DistanceMatrix = zeros(13,13);  %������󣨴��۾���
for k = 2:13                %����ˮվ��һ��ˮվ����
    DistanceMatrix(1,k) = Point_Distance(Point(1),Point(k));
end
for k = 2:13                %һ��ˮվ������ˮվ��Ϊ�����
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

%% ��Ƕ����
function D = Point_Distance(PointS, PointE)
D = sqrt((PointS.x-PointE.x)^2+(PointS.y-PointE.y)^2);
end

function [Tree,Dis]=MinSpanTree(adjMat,sPoint)
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

