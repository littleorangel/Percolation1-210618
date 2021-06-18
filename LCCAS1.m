%% �;���ϵ����·���ӣ�LCCAS----1.0(Clustering Coefficient)
%% �ڽӾ���Ļ�ȡ
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
 %% �ڽӾ���ĳ�������(ȥ����Ϊ0��Ԫ��)
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
[sU, eU] = find(triu(Unbond_matrix) == -1); UnbondLoca = [sU, eU]'; %δ���ӱ߱ߵ�λ�þ���UnbondLoca
Max_NumBond = length(UnbondLoca); %��ʼ�����������ӵ������·��UnbondLoca
%% �����ʼ��ֵR_old
G = graph(A,'upper', 'omitselfloops');
NodeNum = numnodes(G); %�ڵ����
BondNum = numedges(G); %�߸���
DegMean = mean(degree(G)); %MeanDeg�������ƽ����
figure, plot(G,'Layout','force')
% R_old = Percolation(A); %������ֵ
%% ��������·�Ľ׶λ��֣���ȷ��������������ֵ�Ĵ���
Par1 = 1000; %Par1��ÿ�μ�������ֵʱ��·��������
Par2 = 10000; %Par2��������·����
BondNew_ALL = 0; %NewBond�����������������·��
B = A;
%% ��·���ӹ��򣺵��м�������������·��ÿ����Par1����·�����¼������������ֵR_new
R_new = []; %������·�����������ֵ�洢��
T1 = [];
while BondNew_ALL < Par2
    G1 = graph(B,'upper', 'omitselfloops');
    BondNumB = numedges(G1); %�߸���
    BondNewNum  = 0; %������·����
    %% ÿ�ε�������С����ϵ���Ľڵ���ڽӽڵ㹹�ɵ���ͼ������������
    tStart = tic;
    while BondNewNum < Par1 %ȷ������Par1����·��ż���һ������ֵ
        if BondNewNum >= Par1
            break
        end
        B(logical(eye(size(A)))) = 1; %������B�ĶԽ���Ԫ�ػ�Ϊ1
        %% ��֤�ڽ��о���ϵ���ļ���ʱ�������нڵ����С�����ٴ��ڵ���2
        Deg = degree(G1);
        DegMin = min(Deg);
        if DegMin < 2
            DegMinL = find(Deg == DegMin); %�����С�Ƚڵ��λ��DegMinL
            Deg_sort = unique(sort(Deg)); %�ڵ�Ķ�����Deg_sort
            Diff = Deg_sort(2)-Deg_sort(1); %�ڵ�Ķ������������Ĳ���Diff
            B(logical(eye(size(A)))) = 2; %������B�ĶԽ���Ԫ�ػ�Ϊ2
            NodeD = B(DegMinL(1), :); %ȡ��Ŀ��ڵ�DegMinL(1)�ı�����,���ڵ�ı�Ŵ�С��ʧ����ѡȡ
            [sU, eU] = find(NodeD == 0); %��Ŀ��ڵ�δ���ӵĽڵ��λ������[sU, eU]
            Array_Mindeg = degree(G1, eU); %���δ���ӽڵ�Ķȷֲ�����Array_Mindeg
            DegLoca = [Array_Mindeg; eU]'; %DegLoca��δ���ӽڵ�Ķȷֲ�������λ����������һ��Ϊ�ȷֲ��������ڶ���Ϊλ������
            DegLoca_sort = sortrows(DegLoca, 1); %DegLoca_sort��δ���ӽڵ�Ķȷֲ�������λ���������յ�һ�еĶȷֲ���������
            for jj = 1:Diff
                B(DegMinL(1), DegLoca_sort(jj, 2)) = 1;B(DegLoca_sort(jj, 2), DegMinL(1)) = 1; %���յͶȽ�����·������
                G1 = graph(B,'upper', 'omitselfloops');
                BondNewNum = numedges(G1) - BondNumB %���ѭ����������·����
                if BondNewNum >= Par1
                    break
                end
            end
        end
        
        %% ������ڵ�ľ���ϵ��CC������ľ���ϵ��CCNet
        CC = zeros(1, NodeNum); %�����ڵ�ľ���ϵ���洢��
        for i = 1:NodeNum
            N = neighbors(G1, i);
            H = subgraph(G1, N);
            BondH = numedges(H);
            NodeH = numnodes(H);
            CC(i) = 2*BondH/(NodeH*(NodeH-1));
        end
        CC(isnan(CC)) = 0; %������ϵ��������NAN�û�Ϊ0
        CCNet = mean(CC); %��������ľ���ϵ��
        %% ����С����ϵ���Ľڵ���ڽӽڵ㹹�ɵ���ͼ������������
        VecL = [1:1:NodeNum]; %���ڵ��λ������VecL
        CL = [VecL; CC]'; %���ڵ��λ�������;���ϵ��ֵCL
        CCMin = min(CC); CCMinL = find(CC == CCMin); %��ȡ��С����ϵ��ֵCCMin����С����ϵ���Ľڵ�λ��CCMinL
        Neig = neighbors(G1, CCMinL(1)); %��ȡ����С����ϵ�����ڵĽڵ���������Nei
        GH = subgraph(G1, Neig); %��ȡ���ɽڵ�CCMinL(1)���ڽӽڵ�Nei�����ɵ���ͼ
        GH_A = full(adjacency(GH)); %��ȡ����ͼ���ڽӾ���GH_A
        GH_A(logical(eye(size(GH_A)))) = 1;  GH_A = GH_A -1;%����ͼ���ڽӾ���ĶԽ�Ԫ�ػ�Ϊ1����δ���ӱ����ڽӾ����е�ֵ��Ϊ-1
        [sOut, eOut] = find(triu(GH_A) == -1); SubUn = [sOut, eOut]; %��ȡδ���ӱߵ��������SubUn 
        Num1 = round(0.6*size(SubUn,1)); Rand = randperm(size(SubUn, 1), Num1); 
        for i = 1:Num1
            B(Neig(sOut(Rand(i))), Neig(eOut(Rand(i)))) = 1;
            B(Neig(eOut(Rand(i))), Neig(sOut(Rand(i)))) = 1;
            G1 = graph(B,'upper', 'omitselfloops');
            BondNewNum = numedges(G1) - BondNumB %���ѭ����������·����
            if BondNewNum >= Par1
                break
            end
        end
    end
    %% ��·���ӵ�һ��ֵʱ���������������ֵ
    B(logical(eye(size(A)))) = 0; %������B�ĶԽ���Ԫ�ػ�Ϊ0
    BondNew_ALL = numedges(G1) - BondNum %BondNew_ALL������������������·����Bond_num����ʼ����ı߸���
    
    R_new = [R_new, Percolation(B)]; %������ֵ����Percolation�У�����������·���������ֵ
    
    T1 = [T1, toc(tStart)];
end
%% ����������ͼ
G1 = graph(B,'upper', 'omitselfloops');
figure, plot(G1,'Layout','force')
Time = toc; %��ʱ
%% 
PT = [R_new; T1];
pathname = 'D:\desk\�����Խ��\R_new\�������\2\';
filename = '2LCCAS1.mat';
save([pathname, filename], 'PT');

