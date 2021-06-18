clc;clear;tic;
%% %洋葱模型2.0
%%
%提取关系矩阵
data1 = xlsread('C:\Users\11236\Desktop\战略性新兴产业-筛选后\各产业的关系矩阵(2021)\数字创意产业-副本.xlsx','Sheet4');
%提取公司数量
data2 = xlsread('C:\Users\11236\Desktop\战略性新兴产业-筛选后\各产业的关系矩阵(2021)\数字创意产业-副本.xlsx','Sheet5');
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
A = A+A';
A(find(A==2))=[1];
%% 
G = graph(A);
Node_number = numnodes(G); %节点个数
Bond_number = numedges(G); %边个数
figure, plot(G)
R_old = Percolation(A); %计算阈值
B = A;
[start_B, endot_B] = find(triu(B)==1); %交换边后的起点和终点编号
Matri_bond_B = [start_B, endot_B]';  %边的位置矩阵
Success_number = 0; %成功交换的次数的初值值
%% 
for i = 1:1000
    i
    Success_change = 0; %初始化成功交换次数
    %% %成功完成5次边的交换
    C = B; %将交换之前的矩阵先存储一下(B为确定可以交换后的矩阵)
    while Success_change < 5
        Random_array = randperm(Bond_number, 2); %生成随机数组
        Bond1 = [Matri_bond_B(1, Random_array(1)), Matri_bond_B(2, Random_array(1))]; %选取两个待交换的边
        Bond2= [Matri_bond_B(1, Random_array(2)), Matri_bond_B(2, Random_array(2))];
        Bond_VS = [Bond1, Bond2];
        if length(unique(Bond_VS)) == 4  %确保所选的四个节点不同，否则会导致节点的度发生改变
            G_B = graph(B);
            dist1 = distances(G_B, Bond1(1, 1), Bond2(1, 2)); %返回新增加边的最短路径
            dist2 = distances(G_B, Bond2(1, 1), Bond1(1, 2)); %返回新增加边的最短路径
            if dist1 == 1 || dist2 == 1 %确保交换后的边是之前所不存在的，否则会导致节点的度发生改变
            else
                B(Bond1(1, 1), Bond1(1, 2)) = 0; B(Bond1(1, 2), Bond1(1, 1)) = 0; %删除待交换的边
                B(Bond2(1, 1), Bond1(1, 2)) = 0; B(Bond2(1, 2), Bond1(1, 1)) = 0;
                B(Bond1(1, 1), Bond2(1, 2)) = 1; B(Bond2(1, 2), Bond1(1, 1)) = 1; %增加新边
                B(Bond2(1, 1), Bond1(1, 2)) = 1; B(Bond1(1, 2), Bond2(1, 1)) = 1;
                Success_change = Success_change+1;
            end
        else
        end
    end
    %% %成功完成五次交换后计算交换后的新阈值
    R_new = Percolation(B); %带入阈值函数Percolation中，计算交换边后的逾渗阈值
    P_threshold = 0.01; %设置是否交换的阈值
    if R_new > R_old+P_threshold
        [start_B, endot_B] = find(triu(C)==1); %交换边之后更新的起点和终点编号
        Matri_bond_B = [start_B, endot_B]';  %更新边的位置矩阵
        R_old = R_new; %阈值明显提高，成功完成边的交换，新阈值变为下一次迭代的旧阈值
        Success_number = Success_number+1;
    else
        B = C; %如果逾渗阈值增长的不够明显，将矩阵变回未交换边之前
    end
end
toc;