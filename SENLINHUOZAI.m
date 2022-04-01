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
