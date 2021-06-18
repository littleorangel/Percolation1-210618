%% 低合成指标链路增加：LSIAS----1.0
%% 获取邻接矩阵
clc;clear;tic;
data1 = xlsread('C:\Users\11236\Desktop\战略性新兴产业-筛选后\各产业的关系矩阵(2021)-副本\新能源汽车产业-副本.xlsx','Sheet4');%提取关系矩阵
data2 = xlsread('C:\Users\11236\Desktop\战略性新兴产业-筛选后\各产业的关系矩阵(2021)-副本\新能源汽车产业-副本.xlsx','Sheet5');%提取公司数量
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
%% 计算初始阈值R_old
G = graph(A,'upper', 'omitselfloops');
NodeNum = numnodes(G); %节点个数
BondNum = numedges(G); %边个数
MeanDeg = mean(degree(G)); %MeanDeg：网络的平均度
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
        G1 = graph(B,'upper', 'omitselfloops');
        Deg = degree(G1); DegMin = min(Deg); DegMax = max(Deg); %节点的度deg
        Bet = centrality(G1,'betweenness'); BetMin = min(Bet); BetMax = max(Bet); %计算节点的中间中心性betweenness
        Close = centrality(G1,'closeness'); CloseMin = min(Close); CloseMax = max(Close);%网络中节点的中间中心性向量closeness
        CC = zeros(NodeNum, 1);
        for i = 1:NodeNum
            N = neighbors(G1, i);
            H = subgraph(G1, N);
            BondH = numedges(H);
            NodeH = numnodes(H);
            CC(i) = 2*BondH/(NodeH*(NodeH-1));
        end
        CC(isnan(CC)) = 0; CCMin = min(CC); CCMax = max(CC);%计算节点的聚类系数CC
        Si = ((Deg-DegMin)./(DegMax-DegMin) + (Bet-BetMin)./(BetMax-BetMin) + (Close-CloseMin)./(CloseMax-CloseMin) + (CC-CCMin)./(CCMax-CCMin))./4; %合成指标的计算Si
        B(logical(eye(size(B)))) = 1; %将对角线的元素换为1
        
        SiMin = min(Si); %最小合成指标值SiMin
        SiMinL = find(Si == SiMin); %最小合成指标值对应的节点位置SiMinL
        NodeD = B(SiMinL(1), :); %提取目标节点（最小合成指标值）
        [sU, eU] = find(NodeD == 0); %和目标节点未连接的节点的位置向量[sU, eU]
        if isempty(sU)
            c = 2;
            Si_sort = unique(sort(Si));
            SiMin = Si_sort(c);
            SiMinL = find(Si == SiMin);
            NodeD = B(SiMinL(1), :);
            [sU, eU] = find(NodeD == 0); %和目标节点未连接的节点的位置向量[sU, eU]
            if isempty(sU)
                c = c+1;
                continue
            end
        end
        Si_Node = Si(eU); %与目标节点不相连的节点所对应的合成指标值
        SiLo = [Si_Node, eU']; %合成指标-位置矩阵，第一列为合成指标值，第二列为相对应的节点位置
        SiLo_sort = sortrows(SiLo, 1); %将合成指标-位置矩阵按第一列的合成指标值进行排序
        B(SiMinL(1), SiLo_sort(1, 2)) = 1;
        B(SiLo_sort(1, 2), SiMinL(1)) = 1;
        G1 = graph(B,'upper', 'omitselfloops');
        NewBond_num = numedges(G1) - Bond_number_B %完成循环后，新增链路个数
    end
    %% 当新增链路达到一定值时，计算新网络的阈值
    NewBondall = numedges(G1) - BondNum %NewBond：整个新网络新增链路数
    B(logical(eye(size(A)))) = 0; %将矩阵B的对角线元素换为0

    R_new = [R_new, Percolation(B)]; %带入阈值函数Percolation中，计算增加链路后的逾渗阈值

    T1 = [T1, toc(tStart)];
    %% 指定当新网络的最小度大于某值时，则退出整个循环
    if DegMin > (DegMinNew-1)
        break
    end
end

%% 绘制新网络图
G1 = graph(B,'upper', 'omitselfloops');
figure, plot(G1,'Layout','force','EdgeAlpha',0.1)
Time = toc; %计时
%% 
PT = [R_new; T1];
pathname = 'D:\desk\各策略结果\R_new\变量结果\8New\';
filename = '8LSIAS1.mat';
save([pathname, filename], 'PT');



