%% �˶��滮�γ��Ŵ��㷨����ʵ��
% ���ܣ� 
% ������֪�����У�Ѱ�ҵ�������λ������ĵ�
% Created by Yuanlong Zhang
% Date: 2020.10.19

function GeneticAlgorithmSiumlation4

clear;
clc;                        % ����
close all                   % �ر����л�ͼ��
x_max   = 20;               % ��������ֵ��Ĭ�����ú����� x �������� y �Ŀ̶���ͬ
x_min   = 0;                % �������Сֵ
N_point = 10;            	% ������·����������ɸ�����Ҫ�����趨
x_range = repmat([x_min, x_max],[N_point*2,1]);    % N_point���㣬2*N_point�������� repmat��A��[m,n]),����m*n��A
[num_var,~] = size(x_range);
popsize = 200;           	% Ⱥ���С
chromlength = 10;          	% �ַ�������(���峤��)����Ҫ�����������ľ��ȡ�������ȡֵ��Χ�ۺ��ж�
pc = 0.6;               	% ������ʣ�ֻ���������С��pcʱ���Ż�������� һ��ȡ 60~100%
pm = 0.02;                	% ������ʣ�һ�� 0.1~10%
iter_num    = 200;        	% �Ŵ����ĵ�����

% ���ȣ�����ѭ������һ�����еĳ�ʼ·��
isTerminal  = 0;            % ѭ��������־λ����ֵΪ0��1
iter_max = 100;          	% ��ʼ������������Ĭ������ 100��
                            % �����Բ����ҵ����ʵĳ�ʼ�⣬��ʾ��ǰ����·�������޿��н⣬��������С����·��
Counter_i = 1;
while ~isTerminal           % ~ ȡ�����ţ�isTerminal=0������ѭ����1������ѭ��
    % ���룬������������Ʊ���ĳ�ʼȺ��
    pop = initpop(popsize, chromlength, num_var);    
    
    % ����Ŀ��ֵ
    [objvalue, x, y] = calobjvalue(pop, x_range, chromlength);
    
    % ����Ⱥ����ÿ���������Ӧ��
    fitvalue = calfitvalue(x, y, objvalue);
    
    if ~isempty(find(fitvalue)>0) || Counter_i >= iter_max     % ����ȫΪ0����˵���ҵ���һ������ĳ�ֵ��isTerminal =1������ѭ��
         isTerminal =1;                                        % ͬʱ������������Counter_i �����趨ֵʱҲ�������ѭ��       
    else
        Counter_i = Counter_i + 1;
    end    
end

if Counter_i >= iter_max
    error('��ǰ�޿��н⣡��������������ã�')
end

for i = 1:1:iter_num
    % ѡ�񣬸���
    newpop      = selection(pop,fitvalue);     
    
    % ����
    newpop1     = crossover(newpop,pc);                 
    
    % ����
    newpop2     = mutation(newpop1,pm);               
    
    % ����Ŀ�꺯��

    [objvalue, x, y] = calobjvalue(newpop2, x_range, chromlength);
    
    % ����Ⱥ����ÿ���������Ӧ��
    fitvalue = calfitvalue(x, y, objvalue);
    
    % ƽ����Ӧ��
    fitness_ave(i) = calavefitness(fitvalue);         	
    
    % ���Ⱥ������Ӧֵ���ĸ��弰����Ӧֵ
    [bestindividual,bestfit]=best(newpop2,fitvalue); 
    
    % ���ص� best_value ������Ӧ��ֵ�����Ǻ���ֵ
    best_value(i)    = bestfit;        
    
    % ���룬���Ա��������ʮ������
    [x_best(i,:), y_best(i,:)]    = decodechrom(bestindividual, x_range, chromlength);	
    
    % ������Ⱥ
    pop     = newpop2;
end


% ��ͼ
figure(1)
plot(1:iter_num, fitness_ave, 'r', 1:iter_num, best_value, 'b')
grid on
legend('ƽ����Ӧ��', '������Ӧ��')

% ��������ģ�ͻ���
figure(2);
fill([2,4,4,2,2],[0,0,2,2,0],[0,0,0])
hold on
fill([0,2,2,4,4,0,0],[6,6,4,4,8,8,6],[0,0,0])
hold on
fill([2,6,6,2,2],[12,12,14,14,12],[0,0,0])
hold on
fill([0,2,2,0,0],[18,18,20,20,18],[0,0,0])
hold on
fill([6,10,10,6,6],[2,2,6,6,2],[0,0,0])
hold on
fill([8,10,10,8,8],[16,16,18,18,16],[0,0,0])
hold on
fill([10,12,12,10,10],[12,12,14,14,12],[0,0,0])
hold on
fill([12,14,14,12,12],[6,6,8,8,6],[0,0,0])
hold on
fill([12,14,14,12,12],[18,18,20,20,18],[0,0,0])
hold on
fill([14,18,18,16,16,14,14],[12,12,14,14,16,16,12],[0,0,0])
hold on
fill([16,20,20,16,16],[2,2,4,4,2],[0,0,0])
hold on
fill([18,20,20,18,18],[8,8,10,10,8],[0,0,0])
hold on

[z ,index]=max(best_value);                     % �������ֵ����λ��
x_chrom_ = x_best(index,:);                     % ���·�� x ����
y_chrom_ = y_best(index,:);                     % ���·�� y ����
scatter(x_chrom_, y_chrom_, 'bo', 'filled');    % ���·����ɢ��ͼ
plot(x_chrom_,y_chrom_, 'r','Linewidth',2)      % ���ߣ��õ����·��

% ������
ymax = z;
disp(['����Ⱦɫ��Ϊ', num2str([x_chrom_,y_chrom_])])
disp(['���·������Ϊ', num2str(-(ymax-200))])
end


%% 2.1��ʼ��(����)
% initpop.m�����Ĺ�����ʵ��Ⱥ��ĳ�ʼ����popsize��ʾȺ��Ĵ�С��chromlength��ʾȾɫ��ĳ���(��ֵ���ĳ���)��
% rand�������ÿ����ԪΪ {0,1} ����Ϊpopsize������Ϊchromlength�ľ���
% round�Ծ����ÿ����Ԫ����Բ���������롣
% ���ȴ�Сȡ���ڱ����Ķ����Ʊ���ĳ���*����������
% ��ʼ��
function pop=initpop(popsize,chromlength, num_var) 
if num_var == 1                 % ֻ��һ�����������
    pop = round(rand(popsize,chromlength)); 
else                            % ������������
    pop = round(rand(popsize,chromlength*num_var));
end

end

%% 2.2 ����Ŀ�꺯��ֵ
% 2.2.1 calobjvalue.m�����Ĺ�����ʵ��Ŀ�꺯���ļ���
% ʵ��Ŀ�꺯���ļ��㣬�� ��ֵ�� �е���ת��Ϊ ���������
% ����pop ��ʾ��Ⱥ��
% x_range ������ȡֵ��Χ
% chromlength  �����Ʊ���ĳ���
function [objvalue, x, y]=calobjvalue(pop, x_range, chromlength)
[x,y] = decodechrom(pop,x_range,chromlength);
[num_var,~] = size(x_range);
[px,~] = size(pop);
for i = 1:px
    d = 0;
    for j = 1:num_var/2-1
        d = d + sqrt((x(i,j+1)-x(i,j))^2 + (y(i,j+1)-y(i,j))^2);
    end
    objvalue(i,1) = 1 / d;
end
end
%�������ִ����֮��[objvalue , x, y],objvalue
%��һ��200�У�һ�еľ��󣬱�����ÿһ�е�x��y������ֵ,x,y��ÿһ����ʮ����


% 2.2.2 �������Ʊ���ת��Ϊʮ������(2)
% decodechrom.m�����Ĺ����ǽ�Ⱦɫ��(������Ʊ���)ת��Ϊʮ���ƣ�
% ����pop ��ʾ��Ⱥ��
% x_range ������ȡֵ��Χ
% chromlength  �����Ʊ���ĳ���
% (���ڶ���������ԣ�������������������2*chromlengthλ��ʾ��ÿ������chromlengthλ
% ���һ��������1��ʼ����һ��������chromlength+1��ʼ)
function [x_sort, y_sort]=decodechrom(pop,x_range,chromlength)          
[num_var, temp] = size(x_range);    % 
[px,~] = size(pop);
if num_var == 1         % ���ݱ�����ȡֵ��Χ�ж��Ƿ��ж��������num_varΪ1����ʾֻ��һ������
    pop1 = pop;
    x_min = x_range(1);
    x_max = x_range(2);    
    pop2 = x_min + (x_max - x_min) * decodebinary(pop1)/(2^chromlength - 1);  %��popÿ��ת����ʮ����ֵ
else                    % �������ʱ
    pop1 = pop;
    x_min = x_range(1);
    x_max = x_range(21);
    for i = 1:num_var
        popvar = pop1(1:px,(i-1)*chromlength+1:i*chromlength);
        pop2(1:px,i) = x_min + (x_max - x_min) * decodebinary(popvar)/(2^chromlength - 1);
    end
end

% �� pop2 �е����ֱ𸳸������� x �������� y
x = pop2(:,1:num_var/2);
y = pop2(:,num_var/2+1:end);

% Ҫ���ʼ�ӣ�0��0�����ն�Ϊ��20��20��������޸ĳ�ʼ���ն�ֵ
%��ʱx,y�ֱ���һ��200�У�10�еľ���
x(:,1) = 0;
x(:,num_var/2) = 20;
y(:,1) = 0;
y(:,num_var/2) = 20;
%��ʱx,y����ĵ�һ�к����һ�зֱ���0��20��ʮ������һ�б���ʮ����,������һ��·

% ���°���С��������
%�а���С����
x_sort = sort(x,2);     % ���� sort ������������������
y_sort = sort(y,2);

end

% ����������ת��Ϊʮ������
function pop2=decodebinary(pop)
[px,py]=size(pop);                  % ��pop�к�����
for i=1:py
    pop1(:,i)=2.^(py-i).*pop(:,i);  % ���� [2^n 2^(n-1) ... 1] ����������Ȼ�����
end
pop2=sum(pop1,2);                   % ��pop1��ÿ��֮�ͣ��õ���Ӧ��ʮ������

end

%% 2.3 ����������Ӧֵ
% �Ŵ��㷨�ӳ���
% ����������Ӧֵ
function fitvalue = calfitvalue(x, y, objvalue)

%���û���ģ��
%x3��y3��x4��y4 Ϊ�����ϰ���߽�����ڶ�����ɵ��߶Ρ�
x3=[2,4,4,2,2,4,4,0,0,2,2,6,6,2,0,2,2,0,6,10,10,6,8,10,10,8,12,14,14,12,10,12,12,10,12,14,14,12,14,18,18,16,16,14,16,20,20,16,18,20,20,18];
y3=[0,0,2,2,4,4,8,8,6,6,12,12,14,14,18,18,20,20,2,2,6,6,16,16,18,18,6,6,8,8,12,12,14,14,18,18,20,20,12,12,14,14,16,16,2,2,4,4,8,8,10,10];
x4=[4,4,2,2,4,4,0,0,2,2,6,6,2,2,2,2,0,0,10,10,6,6,10,10,8,8,14,14,12,12,12,12,10,10,14,14,12,12,18,18,16,16,14,14,20,20,16,16,20,20,18,18,];
y4=[0,2,2,0,4,8,8,6,6,4,12,14,14,12,18,20,20,18,2,6,6,2,16,18,18,16,6,8,8,6,12,14,14,12,18,20,20,18,12,14,14,16,16,12,2,4,4,2,8,10,10,8];

[row_x,col_x]=size(x);          % ȡ x ���к���
fitvalue = zeros(row_x,1);      % ��ʼ�������� fitvalue
%row_x=200 , col_x = 10;
for j=1:1:row_x
    % ��x1,y1)��ŵ�j�еĵ�1�������ڶ�����
    x1 = x(j,1:col_x-1);
    y1 = y(j,1:col_x-1);
    
    
    % ��x2,y2)��ŵ�j�еĵ�2��������һ����
    x2 = x(j,2:col_x);
    y2 = y(j,2:col_x);
    
    %x1,y1һ��һ��·ʮ�����ǰ�Ÿ���ľ���
    %x2,y2һ��һ��·ʮ����ĺ�Ÿ���ľ���
    % �ж���Ⱥ���ڵ����ɵ��߶��Ƿ��������ϰ���߽��ཻ
    ch = check(x1,y1,x2,y2,x3,y3,x4,y4);   
    
    % ������� [ch] ����Ԫ�صĺͲ�����0,��˵������ĳһ���߶ξ����ϰ���߽硣
    if sum(sum(ch))>0           
        temp = 0;
    else
        temp = objvalue(j);
    end
    fitvalue(j)=temp;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%�������ƣ��߽���Ǽ�⺯�� check.m
%%��ڲ����������߶εĶ��㡣���� �ж��߶�PQ �� XY �Ƿ��ཻ������ P��x1,y1����Q(x2,y2)��X(x3,y3)��Y(x4,y4)
                            %%������������Ϊ��chack(x1,y1,x2,y2,x3,y3,x4,y4)
%%���ڲ��������Ƿ��ཻ 0�����ཻ  1���ཻ
%%˵����
    %%�ú������Լ�������߶��Ƿ��ཻ������ڲ�����������ʱ�򣬿���ͬʱ�������߶��Ƿ��ཻ��
    %%�ú�����ڲ�������Ϊ����������������������Ҳ�������������������뱣֤������ڲ�����ά����ͬ��
    %%��ⷽʽ�����໥����ʵ��ʵ�֡��������߶��໥����ʱ�������߶��ڶ�άƽ�����ཻ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pop]=check(x1,y1,x2,y2,x3,y3,x4,y4)
%x1,y1һ��һ��·ʮ�����ǰ�Ÿ���ľ���
    %x2,y2һ��һ��·ʮ����ĺ�Ÿ���ľ���
[px,py]=size(x1);
p=max(px,py);       % ʵ�������������������Խ��м��
[px,py]=size(x3);
q=max(px,py);
%��ʱpx=1,py=9;p=9;
%px=1,py=52;q=52
%�����໥����ʵ��
for j=1:1:q         % ����߶��ཻ�Ĵ���
    for i=1:1:p     % ���������Ⱥ���ڵ����ɵ��߶��Ƿ��뵱ǰ�ϰ���߽��ཻ
        pabx = x1(i) - x2(i);
        paby = y1(i) - y2(i);
        pacx = x1(i) - x3(j);
        pacy = y1(i) - y3(j);
        padx = x1(i) - x4(j);
        pady = y1(i) - y4(j);
        m = pabx*pacy - pacx*paby;
        n = pabx*pady - padx*paby;
        
        pcdx = x3(j) - x4(j);
        pcdy = y3(j) - y4(j);
        pcax = x3(j) - x1(i);
        pcay = y3(j) - y1(i);
        pcbx = x3(j) - x2(i);
        pcby = y3(j) - y2(i);
        mm = pcdx*pcay - pcax*pcdy;
        nn = pcdx*pcby - pcbx*pcdy;
        
        
        if m*n<=0 && mm*nn<=0     % �ཻ
            pop(j,i)=1;           % ������浽����[pop]�С�
        else                      % ��Ȼ����[pop]��һ��q��p�еľ���
            pop(j,i)=0;
        end
    end
end
end



%% ��ƽ����Ӧ��
function fitness_ave = calavefitness(fitness)
[N ,~] = size(fitness);
fitness_ave = sum(fitness)/N;

end


%% 2.4 ѡ��
% ѡ����
% ѡ����Ʋ����Ǿ�����Щ������Խ�����һ���������в��ö�����ѡ��ѡ�����ַ�������ʵ�֡�
% ���ݷ��� pi=fi/��fi=fi/fsum ��ѡ���裺
% 1�� �ڵ� t �����ɣ�1��ʽ���� fsum �� pi 
% 2�� ���� {0,1} ������� rand( .)���� s=rand( .)*fsum
% 3�� �� ����fi��s ����С�� k ����� k �����屻ѡ��
% 4�� ���� N ��2����3���������õ� N �����壬��Ϊ�� t=t+1 ����Ⱥ
% �Ŵ��㷨�ӳ���
function [newpop]=selection(pop, fitvalue) 
totalfit = sum(fitvalue);                       % ����Ӧֵ֮��
fitvalue = fitvalue/totalfit;                   % �������屻ѡ��ĸ���
fitvalue = cumsum(fitvalue);                    % �� fitvalue=[1 2 3 4]����cumsum(fitvalue)=[1 3 6 10],Ҫ�ۼӣ����̶ķ������ο��Ƿ���ת�õ������� 
[row_p,col_p] = size(pop);                      % row_p*col_p
ms = sort(rand(row_p,1));                       % ��С��������
fitin = 1;
newin = 1;
while newin <= row_p                           	% ѡ��col_p���¸��壬���ظ��������������ܵķ�����̫һ��
    if(ms(newin))<fitvalue(fitin)
        newpop(newin,:) = pop(fitin,:);
        newin = newin + 1;
    else
        fitin = fitin + 1;
    end
end

end


%% 2.5 ����
% ����(crossover)��Ⱥ���е�ÿ������֮�䶼��һ���ĸ��� pc ���棬����������Ӹ����ַ�����ĳһλ��
% ��һ�������ȷ������ʼ���ཻ����������������������еĻ�����������顣���磬����2����������x1��x2Ϊ��
% x1=0100110
% x2=1010001
% ��ÿ������ĵ�3λ��ʼ���棬���ֺ�õ�2���µ��Ӵ�����y1��y2�ֱ�Ϊ��
% y1��0100001
% y2��1010110
% ����2���Ӵ�����ͷֱ������2�����������ĳЩ���������ý��������п����ɸ����������Ӵ���ϳɾ��и����ʺ϶ȵĸ��塣
% ��ʵ�Ͻ������Ŵ��㷨������������ͳ�Ż���������Ҫ�ص�֮һ��
% �Ŵ��㷨�ӳ���
function [newpop]=crossover(pop,pc)                         % pc=0.6
[px,py]=size(pop);
newpop=ones(size(pop));
for i=1:2:px-1                                              % ����Ϊ2���ǽ����ڵ�����������н���
    if(rand<pc)
        cpoint=round(rand*py);
        newpop(i,:)=[pop(i,1:cpoint),pop(i+1,cpoint+1:py)];
        newpop(i+1,:)=[pop(i+1,1:cpoint),pop(i,cpoint+1:py)];
    else
        newpop(i,:)=pop(i,:);
        newpop(i+1,:)=pop(i+1,:);
    end
end

end

%% 2.6 ����
% ����(mutation)�������ͻ���ձ����������Ľ��������С�������ָ�����е�ÿ�������ÿһλ���Ը��� pm ��ת��
% ���ɡ�1����Ϊ��0�������ɡ�0����Ϊ��1����
% �Ŵ��㷨�ı������Կ���ʹ���������������������ܴ��ڵ������ռ䣬��˿��� ��һ���̶��� ���ȫ�����Ž⡣
% �Ŵ��㷨�ӳ���
function [newpop]=mutation(pop,pm)
[px,py]=size(pop);
newpop=ones(size(pop));
for i=1:px
    if(rand<pm)
        mpoint=round(rand*py);              % �����ı������1-10֮��
        if mpoint<=0
            mpoint=1;                   	% ����λ��
        end
        newpop(i,:)=pop(i,:);
        if any(newpop(i,mpoint))==0
            newpop(i,mpoint)=1;
        else
            newpop(i,mpoint)=0;
        end
    else
        newpop(i,:)=pop(i,:);
    end
end

end


%% 2.7 ���Ⱥ����������Ӧֵ�������
% �Ŵ��㷨�ӳ���
% ����� t ��Ⱥ������Ӧֵ����ֵ
function [bestindividual,bestfit]=best(pop,fitvalue)
[px,py]=size(pop);
bestindividual=pop(1,:);
bestfit=fitvalue(1);
for i=2:px
    if fitvalue(i)>bestfit
        bestindividual=pop(i,:);
        bestfit=fitvalue(i);
    end
end

end

