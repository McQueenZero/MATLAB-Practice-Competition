%% ��3�� Ԫ���Զ�����ɭ�ֻ���ģ��
n = 200;     %Ԫ�������С
Plight = 0.000001;  %���缸��
Pgrowth = 0.001;  %��������
UL = [n 1:n-1]; %�Ϻ�������
DR = [2:n 1];   %�º�������
veg = zeros(n,n);        %��ʼ��
% The value of veg:
% empty == 0  
% burning == 1
% green == 2
imh = image(cat(3,veg,veg,veg));    %����ͼ�㲢��ȡͼ����
% veg==?�����߼�ֵ����veg(,)==?�����߼�ֵ
for k = 1:100000
    sum = (veg(UL,:) == 1) + (veg(:,UL) == 1) + (veg(DR,:) == 1) + (veg(:,DR) == 1);
    %���ݹ������ɭ�־����� = �� - �Ż���� + ��������
    veg = 2 * (veg == 2) - ( (veg == 2) & (sum > 0 | (rand(n,n) < Plight)) ) + 2 * ( (veg == 0) & rand(n,n) < Pgrowth);
    imh.CData=cat(3, (veg == 1), (veg == 2), zeros(n)); %���²���������ͼ��
    %                   ��ɫ����    ��ɫ����    ��ɫ����
    drawnow

    %����ѭ����ctrl+C��ֹ
end
