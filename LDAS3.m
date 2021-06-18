%% 低度节点链路增加LDAS----3.0（低度节点的链路增加规则：和低度的未连接节点连接）
%% 获取邻接矩阵A
clc;clear;tic;
data1 = xlsread('C:\Users\11236\Desktop\战略性新兴产业-筛选后\各产业的关系矩阵(2021)-副本\节能环保产业-副本.xlsx','Sheet4');%提取关系矩阵
data2 = xlsread('C:\Users\11236\Desktop\战略性新兴产业-筛选后\各产业的关系矩阵(2021)-副本\节能环保产业-副本.xlsx','Sheet5');%提取公司数量
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
A(find(A==2))=1; %去除权重，也可以称为去除双链接
A(logical(eye(size(A))))=0; %将对角线元素置换为0,排除自链接
G = graph(A,'upper', 'omitselfloops');
Node_num = numnodes(G); %节点个数
Deg = degree(G);
[E, S] = find(Deg==0);
A(E, :) = [];
A(:, E) = [];
%% 
Unbond_matrix = A; Unbond_matrix(logical(eye(size(A)))) = 1; Unbond_matrix = triu(triu(Unbond_matrix)-1); %未连接边的邻接矩阵Unbond_matrix
[sUnBond, eUnBond] = find(triu(Unbond_matrix) == -1); Matri_Unbond = [sUnBond, eUnBond]'; %未连接边边的位置矩阵Matri_Unbond
Max_NumBond = length(Matri_Unbond); %初始网络所能增加的最大链路数Max_NumBond
%% 计算初始阈值R_old
G = graph(A,'upper', 'omitselfloops');
Node_number = numnodes(G); %节点个数
Bond_number = numedges(G); %边个数
figure, plot(G,'Layout','force')
% R_old = Percolation(A); %计算阈值
%% 确定每次增加链路数量（Number_each）
DegMinNew = 50; %DegMinNew：设置增加链路后的新网络所需达到的最低度
Par1 = 1000; %Par1：每次计算新阈值时链路增加数量
Par2 = 10000; %Par2：新增链路总数
NewBondall = 0; %NewBond：整个网络的新增链路数
[start, endot] = findedge(G); Matri_bond = [start, endot]'; %边的位置矩阵Matri_bond
B = A;    
%% 低度节点增加链路以及计算新网络的阈值R_new
R_new = []; %增加链路后新网络的阈值存储器
T1 = [];
while NewBondall < Par2 %计算新阈值的次数
    G1 = graph(B, 'upper', 'omitselfloops');
    Bond_number_B = numedges(G1); %边个数
    NewBond_num  = 0; %新增链路个数
    tStart = tic;
    while NewBond_num < Par1 %确保增加100条链路后才计算一遍新阈值
        Deg = degree(G1); %节点的度deg
        Deg_min = min(Deg); %输出最小度
        loca_MDeg = find(Deg == Deg_min); %输出最小度节点的位置loca_MDeg
        Deg_sort = unique(sort(Deg)); %节点的度排序Deg_sort
        Differ = Deg_sort(2)-Deg_sort(1); %节点的度数倒数两个的差异Differ
        B(logical(eye(size(A)))) = 2; %将矩阵B的对角线元素换为2
        %% 指定当新网络的最小度大于某值时，则退出整个循环
        if Deg_min > (DegMinNew-1)
            break
        end
        %% 将最小度的节点loca_MDeg(j)的度都增加Differ条链路（注：和低度的未连接节点连接）
        Array_Node = B(loca_MDeg(1), :); %取出目标节点loca_MDeg(1)的边向量,按节点的标号大小损失依次选取
        [sUnBond, eUnBond] = find(Array_Node == 0); %和目标节点未连接的节点的位置向量[sUnBond, eUnBond]
        Array_Mindeg = degree(G1, eUnBond); %输出未连接节点的度分布向量Array_Mindeg
        DegLoca = [Array_Mindeg; eUnBond]'; %DegLoca：未连接节点的度分布向量和位置向量，第一列为度分布向量，第二列为位置向量
        DegLoca_sort = sortrows(DegLoca, 1); %DegLoca_sort：未连接节点的度分布向量和位置向量按照第一列的度分布向量排序
        for jj = 1:Differ
            B(loca_MDeg(1), DegLoca_sort(jj, 2)) = 1;B(DegLoca_sort(jj, 2), loca_MDeg(1)) = 1; %按照低度进行链路的增加
            G1 = graph(B,'upper', 'omitselfloops');
            NewBond_num = numedges(G1) - Bond_number_B %完成循环后，新增链路个数
            if NewBond_num >= Par1
                break
            end
        end
    end
    %% 当新增链路达到一定值时，计算新网络的阈值
    B(logical(eye(size(A)))) = 0; %将矩阵B的对角线元素换为0
    NewBondall = numedges(G1) - Bond_number %NewBond：整个新网络新增链路数
    
    R_new = [R_new, Percolation(B)]; %带入阈值函数Percolation中，计算增加链路后的逾渗阈值

    T1 = [T1, toc(tStart)];    
    %% 指定当新网络的最小度大于某值时，则退出整个循环
    if Deg_min > (DegMinNew-1)
        break
    end
end
%% 绘制新网络图
G1 = graph(B,'upper', 'omitselfloops');
figure, plot(G1,'Layout','force','EdgeAlpha',0.1)
Time = toc; %计时
%% 
PT = [R_new; T1];
pathname = 'D:\desk\各策略结果\R_new\变量结果\2\';
filename = '2LDAS3.mat';
save([pathname, filename], 'PT');


%ROld = 0.2940
%result = [0.2272, 0.5231, 0.6525, 0.7178, 0.7661, 0.7918, 0.8035, 0.8189, 0.8435, 0.8636]
%MeanDeg_old = 0.9417; MeanDeg_new = 21.7958