%% ���м���������·����LBAS----1.0
%% ��ȡ�ڽӾ���A
clc;clear;tic;
data1 = xlsread('C:\Users\11236\Desktop\ս�������˲�ҵ-ɸѡ��\����ҵ�Ĺ�ϵ����(2021)-����\���ܻ�����ҵ-����.xlsx','Sheet4');%��ȡ��ϵ����
data2 = xlsread('C:\Users\11236\Desktop\ս�������˲�ҵ-ɸѡ��\����ҵ�Ĺ�ϵ����(2021)-����\���ܻ�����ҵ-����.xlsx','Sheet5');%��ȡ��˾����
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
%% 
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
%% �����ʼ��ֵR_old
G = graph(A,'upper', 'omitselfloops');
Node_num = numnodes(G); %�ڵ����
Bond_num = numedges(G); %�߸���
MeanDeg = mean(degree(G)); %MeanDeg�������ƽ����
figure, plot(G,'Layout','force')
% R_old = Percolation(A); %������ֵ
%% ��������·�Ľ׶λ��֣���ȷ��������������ֵ�Ĵ���
Par1 = 1000; %Par1��ÿ�μ�������ֵʱ��·��������
Par2 = 10000; %Par2��������·����
BondNew_ALL = 0; %NewBond�����������������·��
B = A;
toc;
%% ��·���ӹ��򣺵��м�������������·��ÿ����Par1����·�����¼������������ֵR_new
R_new = []; %������·�����������ֵ�洢��
T1 = [];
while BondNew_ALL < Par2
    G1 = graph(B,'upper', 'omitselfloops');
    BondNumB = numedges(G1); %�߸���
    BondNewNum  = 0; %������·����
    %% ÿ�ε����м���������С�Ľڵ�͵ڶ�С�Ľڵ����ӣ�һ��ֻ����һ����·
    tStart = tic;
    while BondNewNum < Par1 %ȷ������Par1����·��ż���һ������ֵ
        B(logical(eye(size(A)))) = 1; %������B�ĶԽ���Ԫ�ػ�Ϊ1
        BetSco = centrality(G1,'betweenness'); %BetSco�������нڵ���м�����������
        BetLoc = [1:1:Node_num]'; %BetLoc�������нڵ���м�����������Ӧ��λ������
        BetMin = min(BetSco); %BetMin����С�м�������
        BetMinLoc = find(BetSco == BetMin); %BeMinLoc����С�м������Զ�Ӧ��λ������        
        BetLS = [BetLoc, BetSco]; %BetLS���ڵ��м������Ժ����Ӧ��������ġ�������󣬵�һ��Ϊλ���������ڶ���Ϊ�м������Եĵ÷�
        
        NodeD = B(BetMinLoc(1), :); %NodeD��ѡȡĿ��ڵ�BetMinLoc(1)�Ĺ�ϵ����
        [sU, eU] = find(NodeD == 1); %�ҵ�����С�м������Խڵ��Ѿ����ӵĽڵ����sU����������eU
        BetLS(eU, :) = []; %���Ѿ���Ŀ��ڵ������Ľڵ�����ġ�����������ų�
        BetLS_Sort = sortrows(BetLS, 2); %BetLS_Sort�������ġ�������󰴵ڶ��е��м������Խ�������
        
        B(BetMinLoc(1), BetLS_Sort(1, 1)) = 1; B(BetLS_Sort(1, 1), BetMinLoc(1)) = 1; %����С�м������ԵĽڵ�BetMinLoc(1)�͵ڶ�С�Ľڵ�BetLS_Sort(1, 1)����
        G1 = graph(B,'upper', 'omitselfloops');
        BondNewNum = numedges(G1) - BondNumB %���ѭ����������·����
    end
    %% ��·���ӵ�һ��ֵʱ���������������ֵ
    BondNew_ALL = numedges(G1) - Bond_num %BondNew_ALL������������������·����Bond_num����ʼ����ı߸���
    B(logical(eye(size(A)))) = 0; %������B�ĶԽ���Ԫ�ػ�Ϊ0
    
    R_new = [R_new, Percolation(B)]; %������ֵ����Percolation�У�����������·���������ֵ
    
    T1 = [T1, toc(tStart)];
end
%% ����������ͼ
G1 = graph(B,'upper', 'omitselfloops');
figure, plot(G1,'Layout','force','EdgeAlpha',0.1);
Time = toc; %��ʱ
%% 
PT = [R_new; T1];
pathname = 'D:\desk\�����Խ��\R_new\�������\2\';
filename = '2LBAS1.mat';
save([pathname, filename], 'PT');



%ROld = 0.2990
%result = [0.1040, 0.4050, 0.5144, 0.5331, 0.5380, 0.5816, 0.7254, 0.7789, 0.8109, 0.8373]
%MeanDeg_old = 0.9417; MeanDeg_new = 21.7958
