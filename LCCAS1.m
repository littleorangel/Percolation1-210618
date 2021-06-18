%% 低聚类系数链路增加：LCCAS----1.0(Clustering Coefficient)
%% 邻接矩阵的获取
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
 %% 邻接矩阵的初步处理(去除度为0的元素)
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
[sU, eU] = find(triu(Unbond_matrix) == -1); UnbondLoca = [sU, eU]'; %未连接边边的位置矩阵UnbondLoca
Max_NumBond = length(UnbondLoca); %初始网络所能增加的最大链路数UnbondLoca
%% 计算初始阈值R_old
G = graph(A,'upper', 'omitselfloops');
NodeNum = numnodes(G); %节点个数
BondNum = numedges(G); %边个数
DegMean = mean(degree(G)); %MeanDeg：网络的平均度
figure, plot(G,'Layout','force')
% R_old = Percolation(A); %计算阈值
%% 新增总链路的阶段划分，以确定计算新网络阈值的次数
Par1 = 1000; %Par1：每次计算新阈值时链路增加数量
Par2 = 10000; %Par2：新增链路总数
BondNew_ALL = 0; %NewBond：整个网络的新增链路数
B = A;
%% 链路增加规则：低中间中心性增加链路，每增加Par1条链路后重新计算新网络的阈值R_new
R_new = []; %增加链路后新网络的阈值存储器
T1 = [];
while BondNew_ALL < Par2
    G1 = graph(B,'upper', 'omitselfloops');
    BondNumB = numedges(G1); %边个数
    BondNewNum  = 0; %新增链路个数
    %% 每次迭代将最小聚类系数的节点的邻接节点构成的子图的连接数增大
    tStart = tic;
    while BondNewNum < Par1 %确保增加Par1条链路后才计算一遍新阈值
        if BondNewNum >= Par1
            break
        end
        B(logical(eye(size(A)))) = 1; %将矩阵B的对角线元素换为1
        %% 保证在进行聚类系数的计算时，网络中节点的最小度至少大于等于2
        Deg = degree(G1);
        DegMin = min(Deg);
        if DegMin < 2
            DegMinL = find(Deg == DegMin); %输出最小度节点的位置DegMinL
            Deg_sort = unique(sort(Deg)); %节点的度排序Deg_sort
            Diff = Deg_sort(2)-Deg_sort(1); %节点的度数倒数两个的差异Diff
            B(logical(eye(size(A)))) = 2; %将矩阵B的对角线元素换为2
            NodeD = B(DegMinL(1), :); %取出目标节点DegMinL(1)的边向量,按节点的标号大小损失依次选取
            [sU, eU] = find(NodeD == 0); %和目标节点未连接的节点的位置向量[sU, eU]
            Array_Mindeg = degree(G1, eU); %输出未连接节点的度分布向量Array_Mindeg
            DegLoca = [Array_Mindeg; eU]'; %DegLoca：未连接节点的度分布向量和位置向量，第一列为度分布向量，第二列为位置向量
            DegLoca_sort = sortrows(DegLoca, 1); %DegLoca_sort：未连接节点的度分布向量和位置向量按照第一列的度分布向量排序
            for jj = 1:Diff
                B(DegMinL(1), DegLoca_sort(jj, 2)) = 1;B(DegLoca_sort(jj, 2), DegMinL(1)) = 1; %按照低度进行链路的增加
                G1 = graph(B,'upper', 'omitselfloops');
                BondNewNum = numedges(G1) - BondNumB %完成循环后，新增链路个数
                if BondNewNum >= Par1
                    break
                end
            end
        end
        
        %% 计算各节点的聚类系数CC和网络的聚类系数CCNet
        CC = zeros(1, NodeNum); %各个节点的聚类系数存储器
        for i = 1:NodeNum
            N = neighbors(G1, i);
            H = subgraph(G1, N);
            BondH = numedges(H);
            NodeH = numnodes(H);
            CC(i) = 2*BondH/(NodeH*(NodeH-1));
        end
        CC(isnan(CC)) = 0; %将聚类系数向量的NAN置换为0
        CCNet = mean(CC); %计算网络的聚类系数
        %% 将最小聚类系数的节点的邻接节点构成的子图的连接数增大
        VecL = [1:1:NodeNum]; %各节点的位置向量VecL
        CL = [VecL; CC]'; %各节点的位置向量和聚类系数值CL
        CCMin = min(CC); CCMinL = find(CC == CCMin); %获取最小聚类系数值CCMin；最小聚类系数的节点位置CCMinL
        Neig = neighbors(G1, CCMinL(1)); %获取与最小聚类系数相邻的节点坐标向量Nei
        GH = subgraph(G1, Neig); %提取出由节点CCMinL(1)的邻接节点Nei所构成的子图
        GH_A = full(adjacency(GH)); %提取出子图的邻接矩阵GH_A
        GH_A(logical(eye(size(GH_A)))) = 1;  GH_A = GH_A -1;%将子图的邻接矩阵的对角元素换为1，将未连接边在邻接矩阵中的值换为-1
        [sOut, eOut] = find(triu(GH_A) == -1); SubUn = [sOut, eOut]; %获取未连接边的坐标矩阵SubUn 
        Num1 = round(0.6*size(SubUn,1)); Rand = randperm(size(SubUn, 1), Num1); 
        for i = 1:Num1
            B(Neig(sOut(Rand(i))), Neig(eOut(Rand(i)))) = 1;
            B(Neig(eOut(Rand(i))), Neig(sOut(Rand(i)))) = 1;
            G1 = graph(B,'upper', 'omitselfloops');
            BondNewNum = numedges(G1) - BondNumB %完成循环后，新增链路个数
            if BondNewNum >= Par1
                break
            end
        end
    end
    %% 链路增加到一定值时，计算新网络的阈值
    B(logical(eye(size(A)))) = 0; %将矩阵B的对角线元素换为0
    BondNew_ALL = numedges(G1) - BondNum %BondNew_ALL：整个新网络新增链路数；Bond_num：初始网络的边个数
    
    R_new = [R_new, Percolation(B)]; %带入阈值函数Percolation中，计算增加链路后的逾渗阈值
    
    T1 = [T1, toc(tStart)];
end
%% 绘制新网络图
G1 = graph(B,'upper', 'omitselfloops');
figure, plot(G1,'Layout','force')
Time = toc; %计时
%% 
PT = [R_new; T1];
pathname = 'D:\desk\各策略结果\R_new\变量结果\2\';
filename = '2LCCAS1.mat';
save([pathname, filename], 'PT');

