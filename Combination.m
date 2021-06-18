%% 计算组合网络的逾渗阈值
function [PTBN] = Combination(SeriALL, GTot, PT)
for i = 1:9
    for j = 1:9
        if i==j
        else
            S1 = SeriALL{i}; S2 = SeriALL{j};
            CoreNode = intersect(S1(:, 1), S2(:, 1)); %产业共有的核心节点组的公司编码
            CNLab1 = zeros(length(CoreNode), 1); %核心节点组1的网络节点代码
            CNLab2 = zeros(length(CoreNode), 1); %核心节点组2的网络节点代码
            for ii = 1:length(CoreNode)
                CNLab1(ii) = find(S1(:, 1)==CoreNode(ii)); %得到核心节点组
                CNLab2(ii) = find(S2(:, 1)==CoreNode(ii));
            end
            KNL1 = CNLab1;
            KNL2 = CNLab2;
            G1 = GTot{i}; G2 = GTot{j};
            Degree1 = degree(G1, KNL1);
            Degree2 = degree(G2, KNL2);
            B(i, j) = sum(Degree1.*Degree2); %子网络i和j的连接强度系数
        end
    end
end
% 计算子网络的平均连接系数
C = zeros(1, 9);%计算每个产业的连接强度系数
for i = 1:9
    C(i) = sum(B(i, :))/(size(B, 1)-1);
end
CMin = min(C); CMax = max(C);
CSMax = 80000; %连接强度系数的设定最大值
CNor = C./CSMax;
BP_ti = 0;
for i =1:8
    BP_ti = BP_ti+(1/(CNor(i)*PT(i)));
end
PTBN = 1/BP_ti;
end