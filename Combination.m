%% ������������������ֵ
function [PTBN] = Combination(SeriALL, GTot, PT)
for i = 1:9
    for j = 1:9
        if i==j
        else
            S1 = SeriALL{i}; S2 = SeriALL{j};
            CoreNode = intersect(S1(:, 1), S2(:, 1)); %��ҵ���еĺ��Ľڵ���Ĺ�˾����
            CNLab1 = zeros(length(CoreNode), 1); %���Ľڵ���1������ڵ����
            CNLab2 = zeros(length(CoreNode), 1); %���Ľڵ���2������ڵ����
            for ii = 1:length(CoreNode)
                CNLab1(ii) = find(S1(:, 1)==CoreNode(ii)); %�õ����Ľڵ���
                CNLab2(ii) = find(S2(:, 1)==CoreNode(ii));
            end
            KNL1 = CNLab1;
            KNL2 = CNLab2;
            G1 = GTot{i}; G2 = GTot{j};
            Degree1 = degree(G1, KNL1);
            Degree2 = degree(G2, KNL2);
            B(i, j) = sum(Degree1.*Degree2); %������i��j������ǿ��ϵ��
        end
    end
end
% �����������ƽ������ϵ��
C = zeros(1, 9);%����ÿ����ҵ������ǿ��ϵ��
for i = 1:9
    C(i) = sum(B(i, :))/(size(B, 1)-1);
end
CMin = min(C); CMax = max(C);
CSMax = 80000; %����ǿ��ϵ�����趨���ֵ
CNor = C./CSMax;
BP_ti = 0;
for i =1:8
    BP_ti = BP_ti+(1/(CNor(i)*PT(i)));
end
PTBN = 1/BP_ti;
end