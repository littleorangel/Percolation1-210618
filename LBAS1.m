%% 低中间中心性链路增加LBAS----1.0
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
%% 
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
%% 计算初始阈值R_old
G = graph(A,'upper', 'omitselfloops');
Node_num = numnodes(G); %节点个数
Bond_num = numedges(G); %边个数
MeanDeg = mean(degree(G)); %MeanDeg：网络的平均度
figure, plot(G,'Layout','force')
% R_old = Percolation(A); %计算阈值
%% 新增总链路的阶段划分，以确定计算新网络阈值的次数
Par1 = 1000; %Par1：每次计算新阈值时链路增加数量
Par2 = 10000; %Par2：新增链路总数
BondNew_ALL = 0; %NewBond：整个网络的新增链路数
B = A;
toc;
%% 链路增加规则：低中间中心性增加链路，每增加Par1条链路后重新计算新网络的阈值R_new
R_new = []; %增加链路后新网络的阈值存储器
T1 = [];
while BondNew_ALL < Par2
    G1 = graph(B,'upper', 'omitselfloops');
    BondNumB = numedges(G1); %边个数
    BondNewNum  = 0; %新增链路个数
    %% 每次迭代中间中心性最小的节点和第二小的节点连接，一次只增加一条链路
    tStart = tic;
    while BondNewNum < Par1 %确保增加Par1条链路后才计算一遍新阈值
        B(logical(eye(size(A)))) = 1; %将矩阵B的对角线元素换为1
        BetSco = centrality(G1,'betweenness'); %BetSco：网络中节点的中间中心性向量
        BetLoc = [1:1:Node_num]'; %BetLoc：网络中节点的中间中心性所对应的位置向量
        BetMin = min(BetSco); %BetMin：最小中间中心性
        BetMinLoc = find(BetSco == BetMin); %BeMinLoc：最小中间中心性对应的位置向量        
        BetLS = [BetLoc, BetSco]; %BetLS：节点中间中心性和其对应坐标的中心—坐标矩阵，第一列为位置向量，第二列为中间中心性的得分
        
        NodeD = B(BetMinLoc(1), :); %NodeD：选取目标节点BetMinLoc(1)的关系向量
        [sU, eU] = find(NodeD == 1); %找到和最小中间中心性节点已经连接的节点个数sU和坐标向量eU
        BetLS(eU, :) = []; %将已经和目标节点相连的节点从中心—坐标矩阵中排除
        BetLS_Sort = sortrows(BetLS, 2); %BetLS_Sort：将中心—坐标矩阵按第二列的中间中心性进行排序
        
        B(BetMinLoc(1), BetLS_Sort(1, 1)) = 1; B(BetLS_Sort(1, 1), BetMinLoc(1)) = 1; %将最小中间中心性的节点BetMinLoc(1)和第二小的节点BetLS_Sort(1, 1)连接
        G1 = graph(B,'upper', 'omitselfloops');
        BondNewNum = numedges(G1) - BondNumB %完成循环后，新增链路个数
    end
    %% 链路增加到一定值时，计算新网络的阈值
    BondNew_ALL = numedges(G1) - Bond_num %BondNew_ALL：整个新网络新增链路数；Bond_num：初始网络的边个数
    B(logical(eye(size(A)))) = 0; %将矩阵B的对角线元素换为0
    
    R_new = [R_new, Percolation(B)]; %带入阈值函数Percolation中，计算增加链路后的逾渗阈值
    
    T1 = [T1, toc(tStart)];
end
%% 绘制新网络图
G1 = graph(B,'upper', 'omitselfloops');
figure, plot(G1,'Layout','force','EdgeAlpha',0.1);
Time = toc; %计时
%% 
PT = [R_new; T1];
pathname = 'D:\desk\各策略结果\R_new\变量结果\2\';
filename = '2LBAS1.mat';
save([pathname, filename], 'PT');



%ROld = 0.2990
%result = [0.1040, 0.4050, 0.5144, 0.5331, 0.5380, 0.5816, 0.7254, 0.7789, 0.8109, 0.8373]
%MeanDeg_old = 0.9417; MeanDeg_new = 21.7958
