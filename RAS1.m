%% ���������·ģ��1.0
clc;clear;
tic;
%% ��ȡ������ڽӾ���A
data1 = xlsread('C:\Users\11236\Desktop\ս�������˲�ҵ-ɸѡ��\����ҵ�Ĺ�ϵ����(2021)-����\��ط���ҵ-����.xlsx','Sheet4');%��ȡ��ϵ����
data2 = xlsread('C:\Users\11236\Desktop\ս�������˲�ҵ-ɸѡ��\����ҵ�Ĺ�ϵ����(2021)-����\��ط���ҵ-����.xlsx','Sheet5');%��ȡ��˾����
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
A(find(A>1)) = 1; %ȥ��Ȩ�أ�Ҳ���Գ�Ϊȥ��˫����
A(logical(eye(size(A))))=0; %���Խ���Ԫ���û�Ϊ0,�ų�������
G = graph(A,'upper', 'omitselfloops');
Node_num = numnodes(G); %�ڵ����
Deg = degree(G);
[E, S] = find(Deg==0);
A(E, :) = [];
A(:, E) = [];
%% 
Unbond_matrix = A; Unbond_matrix(logical(eye(size(A)))) = 1; Unbond_matrix = triu(triu(Unbond_matrix)-1); %δ���ӱߵ��ڽӾ���Unbond_matrix
[sUnBond, eUnBond] = find(triu(Unbond_matrix) == -1); UnbondLoca = [sUnBond, eUnBond]'; %δ���ӱ߱ߵ�λ�þ���UnbondLoca
Max_NumBond = length(UnbondLoca); %��ʼ�����������ӵ������·��UnbondLoca
%% %�����ʼ��ֵR_old
G = graph(A);
NodeNum = numnodes(G); %�ڵ����
BondNum = numedges(G); %�߸���
% figure, plot(G, 'Layout','force')
% R_old = Percolation(A); %������ֵ
%% ȷ��ÿ��������·������Number_each��
DegMinNew = 22; %DegMinNew������������·�������������ﵽ����Ͷ�
Par1 = 1000; %Par1��ÿ�μ�������ֵʱ��·��������
Par2 = 10000; %Par2��������·����
NewBondall = 0; %NewBond�����������������·��
B = A;   
%% �ͶȽڵ�������·�Լ��������������ֵR_new
R_new = []; %������·�����������ֵ�洢��
T1 = [];
GTot = cell(6, 10);
jjj = 1;
while NewBondall < Par2 %��������ֵ�Ĵ���
    G1 = graph(B, 'upper', 'omitselfloops');
    Bond_number_B = numedges(G1); %�߸���
    NewBond_num  = 0; %������·����
    B(logical(eye(size(B)))) = 1;
    tStart = tic;
    while NewBond_num < Par1 %ȷ������100����·��ż���һ������ֵ
        RandArray1 = randperm(NodeNum, 1);
        NodeD = B(RandArray1, :);
        [sU, eU] = find(NodeD == 0); %��Ŀ��ڵ�δ���ӵĽڵ��λ������[sU, eU]
        if isempty(sU)
            continue
        end
        RandArray2 = randperm(length(sU), 1);
        B(RandArray1, eU(RandArray2)) = 1;
        B(eU(RandArray2), RandArray1) = 1;
        G1 = graph(B,'upper', 'omitselfloops');
        NewBond_num = numedges(G1) - Bond_number_B %���ѭ����������·����
    end
    %% ��������·�ﵽһ��ֵʱ���������������ֵ
    NewBondall = numedges(G1) - BondNum %NewBond������������������·��
    G1 = graph(B,'upper', 'omitselfloops');
    GTot{1, jjj} = G1;
    jjj = jjj+1;
    B(logical(eye(size(A)))) = 0; %������B�ĶԽ���Ԫ�ػ�Ϊ0

%     R_new = [R_new, Percolation(B)]; %������ֵ����Percolation�У�����������·���������ֵ
    
    T1 = [T1, toc(tStart)];
end
% %% ����������ͼ
% G1 = graph(B,'upper', 'omitselfloops');
% figure, plot(G1,'Layout','force','EdgeAlpha',0.1)
% Time = toc; %��ʱ
% 
% %% 
% PT = [R_new; T1];
% pathname = 'D:\desk\�����Խ��\R_new\�������\2\'; %·��
% filename = '2RAS1.mat'; %�ļ���
% save([pathname, filename], 'PT');

