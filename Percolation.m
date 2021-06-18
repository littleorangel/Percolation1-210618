%计算网络的阈值的函数
function [R] = Percolation(A) %输入的矩阵A为对称矩阵，增加阻塞概率后的步骤是在最初的网络上进行的
%% 
A(logical(eye(size(A))))=0; %将对角线元素置换为0,排除自链接
G = graph(A);%去除权重的属性
Node_number = numnodes(G); %节点个数
Bond_number = numedges(G); %边个数
[start, endot] = findedge(G); %起点和终点编号
Matri_bond = [start, endot]'; %边的位置矩阵
%% 
%生成所有最短路径长度矩阵Total_dist
Total_dist = distances(G); %距离存储器
Free_number = length(find(Total_dist>0 & Total_dist<inf))/2; %计算有路径的数量
%求解有路径的边矩阵（边的位置矩阵）
[row_free, col_free] = find(triu(Total_dist)>0 & triu(Total_dist<inf)); %有路径的边的起点和终点
Path_free = [row_free, col_free]'; %自由路径的边的位置矩阵
%% 初始条件的设置
total_percolation = []; %逾渗阈值的存储器
%% 潜在优化，可以选择把阻塞概率划分的更加小，如：p = [0.1:0.05:1]等
p = [0.1:0.1:1]; %阻塞概率值
p_numb = size(p, 2); %增加阻塞概率次数
R = zeros(1, 2); %计算5次最后的逾渗阈值
for R_new_i = 1:2
    for i = 1:1000
        %% %生成随机选择路径
        Random_array = randperm(Free_number, 1); %生成随机数组
        Path_experiment = zeros(2, 1); %生成随机选取的路径坐标的存储器
        Path_experiment(:, 1) = Path_free(:, Random_array); %随机选取的路径坐标
        %取一条路径计算逾渗阈值
        ii = Path_experiment(1, 1); %起始节点
        jj = Path_experiment(2, 1); %终端节点
        %%
        %循环求解自由路径数
        n_h = 0;  %自由路径次数初始值
        Total_number = 0;
        Blockbond_mark = []; %标记被阻塞链路的存储器
        for j = 1:p_numb
            B = A; %当概率由P1变为P2时，保证增加阻塞概率后，在原始矩阵的基础上变化
            Bond_number_dele = round(p(j)*Bond_number); %删除边数量
            Name_random = randperm(Bond_number, Bond_number_dele); %生成随机数组
            Blockbond_mark = unique([Blockbond_mark, Name_random]); %将已经选取被阻塞的边标记存储
            for j=1:Bond_number_dele
                B(Matri_bond(1, Name_random(j)), Matri_bond(2, Name_random(j))) = 0; %在阻塞概率为p时，删除边
                B(Matri_bond(2, Name_random(j)), Matri_bond(1, Name_random(j))) = 0;
            end
            G1 = graph(B);
            dist1 = distances(G1, ii, jj); %返回节点1到节点5的最短路径
            if dist1 == inf
                n_h = n_h+0;
                Total_number = Total_number+1;
            else
                n_h = n_h+1;
                Total_number = Total_number+1;
            end
            Mark_number = length(Blockbond_mark);
            if Mark_number == Bond_number
                break
            end
        end
        percolation = n_h/Total_number;
        total_percolation = [total_percolation, percolation];
    end
    %% %对所有实验获得的阈值求平均值获得网络的逾渗阈值
    R_new = mean(total_percolation)
    R(R_new_i) = R_new;
    total_percolation = [];
end
R = mean(R); 
end
