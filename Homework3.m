%%-------------------------------------------------------------------------
% ���ߣ�   ������
% ���ڣ�   2021��4��7��
% ˵����   ��ƽ����100�����͹����αհ����ö�������
% �Ա���룬����PO�������㷨�ο�https://blog.csdn.net/viafcccy/article/details/87483567
%%-------------------------------------------------------------------------
%% ��ƽ����100�����͹����αհ����ö�������
clc, clear, close all

n = 100;
% rng(3);
% 200*200������̬�ֲ�
x = -100 + (100 - (-100))*randn(n,1);
y = -100 + (100 - (-100))*randn(n,1);
% 200*200��������ֲ�
% x = -100 + (100 - (-100))*rand(n,1);
% y = -100 + (100 - (-100))*rand(n,1);

% ---------MATLAB�Դ�����ʵ��
% dt=delaunayTriangulation(x,y);
% k = convexHull(dt);
% plot(x,y,'*');
% hold on;
% plot(x(k),y(k),'r');

% ---------�Ա���루Grahamɨ�跨��
% ��ṹ�帳ֵ
for m = 1:n
    Point(m).x = x(m);
    Point(m).y = y(m);
end

[ysort,ysortindex] = sort(y);  %����������,����������С����Ϊ��ʼ��
xmove = x-x(ysortindex(1));
ymove = y-ysort(1); %����������С����Ϊԭ��
subplot(2,2,1)
plot(xmove,ymove,'*');   %ƽ�ƺ�ĵ�ͼ
title('ƽ�ƺ�ĵ�ͼ')

for m = 1:n
    PointMove(m).x = xmove(m);
    PointMove(m).y = ymove(m);
    if ymove(m) == 0
        y0index = m;    %ԭ������
    end
    PAarray(m) = Point_PolarAngle(PointMove(m));  %��������
end
PAarray(y0index) = -1; %ԭ��ķ�����Ϊ-1�㣬��Ϊԭ��ı�־���������Ǿ�Ϊ�Ǹ���
[PAsort, PAindex] =sort(PAarray);   %������С��������

subplot(2,2,2)
plot(x,y,'*'),hold on
plot(x(PAindex),y(PAindex),'r');     %����ʱ�뷽������
title('����ʱ�뷽������')

% �����ѭ������˰�Χ���е�ıջ�����
for k = 3:n-1
    SPAi = size(PAindex);
    if k+1 > SPAi(2)
        break;
    end
    A = Point_Angle(Point(PAindex(k-1)),Point(PAindex(k)),Point(PAindex(k+1)));
    if A == 1
        PAindex(k)=[];  %��ת˳ʱ����������
        k = k - 1;  
    end 
end
PAindexCloseLine = [PAindex y0index];

subplot(2,2,3)
plot(x,y,'*'),hold on
plot(x(PAindexCloseLine),y(PAindexCloseLine),'r')     %�ջ�����
title('�ջ���������������')

% �����ѭ�����͹��
k = 2;
while 1
    SPAi = size(PAindex);
    if k+1 > SPAi(2)
        for m = 2:SJde(2)-1
            %��ѯÿ�ֽ�������һ���ڰ����Ҷ�Ӧ�ж�ֵ����0
            %�˴������ж�ֵ��Ϊ��У��
            J(m) = Point_Angle(Point(PAindex(m-1)),Point(PAindex(m)),Point(PAindex(m+1)));
            if J(m)==1 %��1�������ұ�/˳ʱ�룩���������
                break;
            end
        end
        k = 2;
        if J == 0 %ȫ0������while loop
            break;
        end
    end
    A = Point_Angle(Point(PAindex(k-1)),Point(PAindex(k)),Point(PAindex(k+1)));
    J(k) = Point_Angle(Point(PAindex(k-1)),Point(PAindex(k)),Point(PAindex(k+1)));
    if A == 1
        PAindex(k)=[];  %��ת˳ʱ���������� 
        J(k)=[];
        k = k - 1;  
    end
    SJde = size(J);
    k = k + 1;
end

PAindex = [PAindex y0index];    %�հ�

subplot(2,2,4)
plot(x,y,'*'),hold on
plot(x(PAindex),y(PAindex),'r')     %͹������
title('͹����αհ�')

TheArea = polyarea(x(PAindex),y(PAindex));
% disp(['͹����ε������',num2str(TheArea)])
figure
plot(x,y,'*'),hold on
plot(x(PAindex),y(PAindex),'r')    %͹������
title(['͹����αհ������Ϊ',num2str(TheArea)])


%% ��Ƕ����
function D = Point_Distance(PointS, PointE)
D = sqrt((PointS.x-PointE.x)^2+(PointS.y-PointE.y)^2);
end

function PA = Point_PolarAngle(PointA)   %�����������
PA = rad2deg(atan(PointA.y/PointA.x));
if PA < 0
    PA = 180 + PA;
end
end

function Bool = Point_Angle(PointS, PointC, PointE)    %�ж�ֱ������
Sidecur = [PointC.x PointC.y]-[PointS.x PointS.y];  %��ǰ�ߵ�����
Sideaft = [PointE.x PointE.y]-[PointC.x PointC.y];  %��һ���ߵ�����
Crorslt = Sidecur(1)*Sideaft(2)-Sidecur(2)*Sideaft(1);
if Crorslt > 0  %���ݲ�������ж���/�ң���/˳ʱ�룩
    Bool = 0;   %��ת
else
    Bool = 1;   %��ת
end
end
