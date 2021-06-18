%�����������ֵ�ĺ���
function [R] = Percolation(A) %����ľ���AΪ�Գƾ��������������ʺ�Ĳ�����������������Ͻ��е�
%% 
A(logical(eye(size(A))))=0; %���Խ���Ԫ���û�Ϊ0,�ų�������
G = graph(A);%ȥ��Ȩ�ص�����
Node_number = numnodes(G); %�ڵ����
Bond_number = numedges(G); %�߸���
[start, endot] = findedge(G); %�����յ���
Matri_bond = [start, endot]'; %�ߵ�λ�þ���
%% 
%�����������·�����Ⱦ���Total_dist
Total_dist = distances(G); %����洢��
Free_number = length(find(Total_dist>0 & Total_dist<inf))/2; %������·��������
%�����·���ı߾��󣨱ߵ�λ�þ���
[row_free, col_free] = find(triu(Total_dist)>0 & triu(Total_dist<inf)); %��·���ıߵ������յ�
Path_free = [row_free, col_free]'; %����·���ıߵ�λ�þ���
%% ��ʼ����������
total_percolation = []; %������ֵ�Ĵ洢��
%% Ǳ���Ż�������ѡ����������ʻ��ֵĸ���С���磺p = [0.1:0.05:1]��
p = [0.1:0.1:1]; %��������ֵ
p_numb = size(p, 2); %�����������ʴ���
R = zeros(1, 2); %����5������������ֵ
for R_new_i = 1:2
    for i = 1:1000
        %% %�������ѡ��·��
        Random_array = randperm(Free_number, 1); %�����������
        Path_experiment = zeros(2, 1); %�������ѡȡ��·������Ĵ洢��
        Path_experiment(:, 1) = Path_free(:, Random_array); %���ѡȡ��·������
        %ȡһ��·������������ֵ
        ii = Path_experiment(1, 1); %��ʼ�ڵ�
        jj = Path_experiment(2, 1); %�ն˽ڵ�
        %%
        %ѭ���������·����
        n_h = 0;  %����·��������ʼֵ
        Total_number = 0;
        Blockbond_mark = []; %��Ǳ�������·�Ĵ洢��
        for j = 1:p_numb
            B = A; %��������P1��ΪP2ʱ����֤�����������ʺ���ԭʼ����Ļ����ϱ仯
            Bond_number_dele = round(p(j)*Bond_number); %ɾ��������
            Name_random = randperm(Bond_number, Bond_number_dele); %�����������
            Blockbond_mark = unique([Blockbond_mark, Name_random]); %���Ѿ�ѡȡ�������ı߱�Ǵ洢
            for j=1:Bond_number_dele
                B(Matri_bond(1, Name_random(j)), Matri_bond(2, Name_random(j))) = 0; %����������Ϊpʱ��ɾ����
                B(Matri_bond(2, Name_random(j)), Matri_bond(1, Name_random(j))) = 0;
            end
            G1 = graph(B);
            dist1 = distances(G1, ii, jj); %���ؽڵ�1���ڵ�5�����·��
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
    %% %������ʵ���õ���ֵ��ƽ��ֵ��������������ֵ
    R_new = mean(total_percolation)
    R(R_new_i) = R_new;
    total_percolation = [];
end
R = mean(R); 
end
