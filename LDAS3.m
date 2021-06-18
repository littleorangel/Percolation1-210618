%% �ͶȽڵ���·����LDAS----3.0���ͶȽڵ����·���ӹ��򣺺͵Ͷȵ�δ���ӽڵ����ӣ�
%% ��ȡ�ڽӾ���A
clc;clear;tic;
data1 = xlsread('C:\Users\11236\Desktop\ս�������˲�ҵ-ɸѡ��\����ҵ�Ĺ�ϵ����(2021)-����\���ܻ�����ҵ-����.xlsx','Sheet4');%��ȡ��ϵ����
data2 = xlsread('C:\Users\11236\Desktop\ս�������˲�ҵ-ɸѡ��\����ҵ�Ĺ�ϵ����(2021)-����\���ܻ�����ҵ-����.xlsx','Sheet5');%��ȡ��˾����
compare_numb = size(data2,1);%��˾����
A = zeros(compare_numb, compare_numb); %�����ڽӾ�����ܾ���
row1 = size(data1,1);
for i = 1:row1 %�ڲ�ѭ���ǶԵ���һ�н���ѭ����iָ����������ÿһ�н���ѭ��
    a1 = data1(i, :);%��ȡdata1�����е�i�е����ݣ����У�ָ��i�е�������
    a1(a1==0) = [];%���������е�0ֵȥ����������Ҫ������
    lie = size(a1, 2);%bָa1������飨���������ݵģ�������
    q = lie-1;%��Ϊ��ѭ��������Ҫ��ԭ����������1��������6�У��������1�������֣�1,7�������
    for k = 1:q
        w = k-1;
        for z = 1:lie-1%��forѭ���ǹ̶���1�У����㣨1,1������1,6��
            hzb = a1(1, 1+w);%wҪ��ȡ0����֤��ȡ����1,1��
            zzb = a1(1, z+w+1);%��֤�ӣ�1,1������1,2���ڵ���1,6��
            A(hzb, zzb) = 1;
        end
        lie = lie-1;
    end
end
%% �ڽӾ���ĳ�������
A = A+A';
A(find(A==2))=1; %ȥ��Ȩ�أ�Ҳ���Գ�Ϊȥ��˫����
A(logical(eye(size(A))))=0; %���Խ���Ԫ���û�Ϊ0,�ų�������
G = graph(A,'upper', 'omitselfloops');
Node_num = numnodes(G); %�ڵ����
Deg = degree(G);
[E, S] = find(Deg==0);
A(E, :) = [];
A(:, E) = [];
%% 
Unbond_matrix = A; Unbond_matrix(logical(eye(size(A)))) = 1; Unbond_matrix = triu(triu(Unbond_matrix)-1); %δ���ӱߵ��ڽӾ���Unbond_matrix
[sUnBond, eUnBond] = find(triu(Unbond_matrix) == -1); Matri_Unbond = [sUnBond, eUnBond]'; %δ���ӱ߱ߵ�λ�þ���Matri_Unbond
Max_NumBond = length(Matri_Unbond); %��ʼ�����������ӵ������·��Max_NumBond
%% �����ʼ��ֵR_old
G = graph(A,'upper', 'omitselfloops');
Node_number = numnodes(G); %�ڵ����
Bond_number = numedges(G); %�߸���
figure, plot(G,'Layout','force')
% R_old = Percolation(A); %������ֵ
%% ȷ��ÿ��������·������Number_each��
DegMinNew = 50; %DegMinNew������������·�������������ﵽ����Ͷ�
Par1 = 1000; %Par1��ÿ�μ�������ֵʱ��·��������
Par2 = 10000; %Par2��������·����
NewBondall = 0; %NewBond�����������������·��
[start, endot] = findedge(G); Matri_bond = [start, endot]'; %�ߵ�λ�þ���Matri_bond
B = A;    
%% �ͶȽڵ�������·�Լ��������������ֵR_new
R_new = []; %������·�����������ֵ�洢��
T1 = [];
while NewBondall < Par2 %��������ֵ�Ĵ���
    G1 = graph(B, 'upper', 'omitselfloops');
    Bond_number_B = numedges(G1); %�߸���
    NewBond_num  = 0; %������·����
    tStart = tic;
    while NewBond_num < Par1 %ȷ������100����·��ż���һ������ֵ
        Deg = degree(G1); %�ڵ�Ķ�deg
        Deg_min = min(Deg); %�����С��
        loca_MDeg = find(Deg == Deg_min); %�����С�Ƚڵ��λ��loca_MDeg
        Deg_sort = unique(sort(Deg)); %�ڵ�Ķ�����Deg_sort
        Differ = Deg_sort(2)-Deg_sort(1); %�ڵ�Ķ������������Ĳ���Differ
        B(logical(eye(size(A)))) = 2; %������B�ĶԽ���Ԫ�ػ�Ϊ2
        %% ָ�������������С�ȴ���ĳֵʱ�����˳�����ѭ��
        if Deg_min > (DegMinNew-1)
            break
        end
        %% ����С�ȵĽڵ�loca_MDeg(j)�Ķȶ�����Differ����·��ע���͵Ͷȵ�δ���ӽڵ����ӣ�
        Array_Node = B(loca_MDeg(1), :); %ȡ��Ŀ��ڵ�loca_MDeg(1)�ı�����,���ڵ�ı�Ŵ�С��ʧ����ѡȡ
        [sUnBond, eUnBond] = find(Array_Node == 0); %��Ŀ��ڵ�δ���ӵĽڵ��λ������[sUnBond, eUnBond]
        Array_Mindeg = degree(G1, eUnBond); %���δ���ӽڵ�Ķȷֲ�����Array_Mindeg
        DegLoca = [Array_Mindeg; eUnBond]'; %DegLoca��δ���ӽڵ�Ķȷֲ�������λ����������һ��Ϊ�ȷֲ��������ڶ���Ϊλ������
        DegLoca_sort = sortrows(DegLoca, 1); %DegLoca_sort��δ���ӽڵ�Ķȷֲ�������λ���������յ�һ�еĶȷֲ���������
        for jj = 1:Differ
            B(loca_MDeg(1), DegLoca_sort(jj, 2)) = 1;B(DegLoca_sort(jj, 2), loca_MDeg(1)) = 1; %���յͶȽ�����·������
            G1 = graph(B,'upper', 'omitselfloops');
            NewBond_num = numedges(G1) - Bond_number_B %���ѭ����������·����
            if NewBond_num >= Par1
                break
            end
        end
    end
    %% ��������·�ﵽһ��ֵʱ���������������ֵ
    B(logical(eye(size(A)))) = 0; %������B�ĶԽ���Ԫ�ػ�Ϊ0
    NewBondall = numedges(G1) - Bond_number %NewBond������������������·��
    
    R_new = [R_new, Percolation(B)]; %������ֵ����Percolation�У�����������·���������ֵ

    T1 = [T1, toc(tStart)];    
    %% ָ�������������С�ȴ���ĳֵʱ�����˳�����ѭ��
    if Deg_min > (DegMinNew-1)
        break
    end
end
%% ����������ͼ
G1 = graph(B,'upper', 'omitselfloops');
figure, plot(G1,'Layout','force','EdgeAlpha',0.1)
Time = toc; %��ʱ
%% 
PT = [R_new; T1];
pathname = 'D:\desk\�����Խ��\R_new\�������\2\';
filename = '2LDAS3.mat';
save([pathname, filename], 'PT');


%ROld = 0.2940
%result = [0.2272, 0.5231, 0.6525, 0.7178, 0.7661, 0.7918, 0.8035, 0.8189, 0.8435, 0.8636]
%MeanDeg_old = 0.9417; MeanDeg_new = 21.7958