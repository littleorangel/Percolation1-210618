%% 随机增加链路模型1.0
clc;clear;
tic;
%% 获取网络的邻接矩阵A
data1 = xlsread('C:\Users\11236\Desktop\战略性新兴产业-筛选后\各产业的关系矩阵(2021)-副本\相关服务业-副本.xlsx','Sheet4');%提取关系矩阵
data2 = xlsread('C:\Users\11236\Desktop\战略性新兴产业-筛选后\各产业的关系矩阵(2021)-副本\相关服务业-副本.xlsx','Sheet5');%提取公司数量
compare_numb = size(data2,1);%公司数量
A = zeros(compare_numb, compare_numb); %建立邻接矩阵接受矩阵
row1 = size(data1,1);
for i = 1:row1 %内部循环是对单个一行进行循环，i指的行数，对每一行进行循环
    a1 = data1(i, :);%提取data1数据中第i行的内容，其中：指第i行的所有列
    a1(a1==0) = [];%将此数组中的0值去掉，留下需要的数据
    lie = size(a1, 2);%b指a1这个数组（都是有内容的）的列数
    q = lie-1;%因为被循环的列数要比原本的列数少1，比如有6列，如果不减1，则会出现（1,7）的情况
    for k = 1:q
        w = k-1;
        for z = 1:lie-1%此for循环是固定第1行，先算（1,1）到（1,6）
            hzb = a1(1, 1+w);%w要先取0，保证能取到（1,1）
            zzb = a1(1, z+w+1);%保证从（1,1）到（1,2）在到（1,6）
            A(hzb, zzb) = 1;
        end
        lie = lie-1;
    end
end
%% 邻接矩阵的初步处理
A = A+A';
A(find(A>1)) = 1; %去除权重，也可以称为去除双链接
A(logical(eye(size(A))))=0; %将对角线元素置换为0,排除自链接
G = graph(A,'upper', 'omitselfloops');
Node_num = numnodes(G); %节点个数
Deg = degree(G);
[E, S] = find(Deg==0);
A(E, :) = [];
A(:, E) = [];
%% 
Unbond_matrix = A; Unbond_matrix(logical(eye(size(A)))) = 1; Unbond_matrix = triu(triu(Unbond_matrix)-1); %未连接边的邻接矩阵Unbond_matrix
[sUnBond, eUnBond] = find(triu(Unbond_matrix) == -1); UnbondLoca = [sUnBond, eUnBond]'; %未连接边边的位置矩阵UnbondLoca
Max_NumBond = length(UnbondLoca); %初始网络所能增加的最大链路数UnbondLoca
%% %计算初始阈值R_old
G = graph(A);
NodeNum = numnodes(G); %节点个数
BondNum = numedges(G); %边个数
% figure, plot(G, 'Layout','force')
% R_old = Percolation(A); %计算阈值
%% 确定每次增加链路数量（Number_each）
DegMinNew = 22; %DegMinNew：设置增加链路后的新网络所需达到的最低度
Par1 = 1000; %Par1：每次计算新阈值时链路增加数量
Par2 = 10000; %Par2：新增链路总数
NewBondall = 0; %NewBond：整个网络的新增链路数
B = A;   
%% 低度节点增加链路以及计算新网络的阈值R_new
R_new = []; %增加链路后新网络的阈值存储器
T1 = [];
GTot = cell(6, 10);
jjj = 1;
while NewBondall < Par2 %计算新阈值的次数
    G1 = graph(B, 'upper', 'omitselfloops');
    Bond_number_B = numedges(G1); %边个数
    NewBond_num  = 0; %新增链路个数
    B(logical(eye(size(B)))) = 1;
    tStart = tic;
    while NewBond_num < Par1 %确保增加100条链路后才计算一遍新阈值
        RandArray1 = randperm(NodeNum, 1);
        NodeD = B(RandArray1, :);
        [sU, eU] = find(NodeD == 0); %和目标节点未连接的节点的位置向量[sU, eU]
        if isempty(sU)
            continue
        end
        RandArray2 = randperm(length(sU), 1);
        B(RandArray1, eU(RandArray2)) = 1;
        B(eU(RandArray2), RandArray1) = 1;
        G1 = graph(B,'upper', 'omitselfloops');
        NewBond_num = numedges(G1) - Bond_number_B %完成循环后，新增链路个数
    end
    %% 当新增链路达到一定值时，计算新网络的阈值
    NewBondall = numedges(G1) - BondNum %NewBond：整个新网络新增链路数
    G1 = graph(B,'upper', 'omitselfloops');
    GTot{1, jjj} = G1;
    jjj = jjj+1;
    B(logical(eye(size(A)))) = 0; %将矩阵B的对角线元素换为0

%     R_new = [R_new, Percolation(B)]; %带入阈值函数Percolation中，计算增加链路后的逾渗阈值
    
    T1 = [T1, toc(tStart)];
end
% %% 绘制新网络图
% G1 = graph(B,'upper', 'omitselfloops');
% figure, plot(G1,'Layout','force','EdgeAlpha',0.1)
% Time = toc; %计时
% 
% %% 
% PT = [R_new; T1];
% pathname = 'D:\desk\各策略结果\R_new\变量结果\2\'; %路径
% filename = '2RAS1.mat'; %文件名
% save([pathname, filename], 'PT');

