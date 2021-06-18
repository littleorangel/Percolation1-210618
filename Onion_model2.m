clc;clear;tic;
%% %���ģ��2.0
%%
%��ȡ��ϵ����
data1 = xlsread('C:\Users\11236\Desktop\ս�������˲�ҵ-ɸѡ��\����ҵ�Ĺ�ϵ����(2021)\���ִ����ҵ-����.xlsx','Sheet4');
%��ȡ��˾����
data2 = xlsread('C:\Users\11236\Desktop\ս�������˲�ҵ-ɸѡ��\����ҵ�Ĺ�ϵ����(2021)\���ִ����ҵ-����.xlsx','Sheet5');
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
A = A+A';
A(find(A==2))=[1];
%% 
G = graph(A);
Node_number = numnodes(G); %�ڵ����
Bond_number = numedges(G); %�߸���
figure, plot(G)
R_old = Percolation(A); %������ֵ
B = A;
[start_B, endot_B] = find(triu(B)==1); %�����ߺ�������յ���
Matri_bond_B = [start_B, endot_B]';  %�ߵ�λ�þ���
Success_number = 0; %�ɹ������Ĵ����ĳ�ֵֵ
%% 
for i = 1:1000
    i
    Success_change = 0; %��ʼ���ɹ���������
    %% %�ɹ����5�αߵĽ���
    C = B; %������֮ǰ�ľ����ȴ洢һ��(BΪȷ�����Խ�����ľ���)
    while Success_change < 5
        Random_array = randperm(Bond_number, 2); %�����������
        Bond1 = [Matri_bond_B(1, Random_array(1)), Matri_bond_B(2, Random_array(1))]; %ѡȡ�����������ı�
        Bond2= [Matri_bond_B(1, Random_array(2)), Matri_bond_B(2, Random_array(2))];
        Bond_VS = [Bond1, Bond2];
        if length(unique(Bond_VS)) == 4  %ȷ����ѡ���ĸ��ڵ㲻ͬ������ᵼ�½ڵ�Ķȷ����ı�
            G_B = graph(B);
            dist1 = distances(G_B, Bond1(1, 1), Bond2(1, 2)); %���������ӱߵ����·��
            dist2 = distances(G_B, Bond2(1, 1), Bond1(1, 2)); %���������ӱߵ����·��
            if dist1 == 1 || dist2 == 1 %ȷ��������ı���֮ǰ�������ڵģ�����ᵼ�½ڵ�Ķȷ����ı�
            else
                B(Bond1(1, 1), Bond1(1, 2)) = 0; B(Bond1(1, 2), Bond1(1, 1)) = 0; %ɾ���������ı�
                B(Bond2(1, 1), Bond1(1, 2)) = 0; B(Bond2(1, 2), Bond1(1, 1)) = 0;
                B(Bond1(1, 1), Bond2(1, 2)) = 1; B(Bond2(1, 2), Bond1(1, 1)) = 1; %�����±�
                B(Bond2(1, 1), Bond1(1, 2)) = 1; B(Bond1(1, 2), Bond2(1, 1)) = 1;
                Success_change = Success_change+1;
            end
        else
        end
    end
    %% %�ɹ������ν�������㽻���������ֵ
    R_new = Percolation(B); %������ֵ����Percolation�У����㽻���ߺ��������ֵ
    P_threshold = 0.01; %�����Ƿ񽻻�����ֵ
    if R_new > R_old+P_threshold
        [start_B, endot_B] = find(triu(C)==1); %������֮����µ������յ���
        Matri_bond_B = [start_B, endot_B]';  %���±ߵ�λ�þ���
        R_old = R_new; %��ֵ������ߣ��ɹ���ɱߵĽ���������ֵ��Ϊ��һ�ε����ľ���ֵ
        Success_number = Success_number+1;
    else
        B = C; %���������ֵ�����Ĳ������ԣ���������δ������֮ǰ
    end
end
toc;