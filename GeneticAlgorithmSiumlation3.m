%% �˶��滮�γ��Ŵ��㷨����ʵ��
% ���ܣ� 
% ������֪�����У�Ѱ�ҵ�������λ������ĵ�
% Created by Yuanlong Zhang
% Date: 2020.10.19

function GeneticAlgorithmSiumlation3

clear;
clc;                        % ����
close all                   % �ر����л�ͼ����
x_max   = 10;               % �Ա���x �����ֵ
x_min   = -10;              % �Ա���x ����Сֵ
x_range = [x_min, x_max;x_min, x_max];
[num_var,~] = size(x_range);
popsize = 200;           	% Ⱥ���С
chromlength = 10;          	% �ַ�������(���峤��)����Ҫ�����������ľ��ȡ�������ȡֵ��Χ�ۺ��ж�
pc = 0.6;               	% ������ʣ�ֻ���������С��pcʱ���Ż�������� һ��ȡ 60~100%
pm = 0.02;                	% ������ʣ�һ�� 0.1~10%
iter_num    = 1000;         % �Ŵ�����

pop = initpop(popsize, chromlength, num_var);             % ������������Ʊ���ĳ�ʼȺ��
for i = 1:1:iter_num
    % ����Ŀ�꺯��
    objvalue    = calobjvalue(pop, x_range, chromlength);      
    
    % ����Ⱥ����ÿ���������Ӧ��
    fitvalue    = calfitvalue(objvalue);          
    
    % ѡ�񣬸���
    newpop      = selection(pop,fitvalue);    
    
    % ����
    newpop1     = crossover(newpop,pc);            
    
    % ����
    newpop2     = mutation(newpop1,pm);                 
    
    % ����Ŀ�꺯��
    objvalue    = calobjvalue(newpop2, x_range, chromlength); 
    
    % ����Ⱥ����ÿ���������Ӧ��
    fitvalue    = calfitvalue(objvalue);            
    
    % ƽ����Ӧ��
    fitness_ave(i) = calavefitness(fitvalue);         	
    
    % ���Ⱥ������Ӧֵ���ĸ��弰����Ӧֵ
    [bestindividual,bestfit]=best(newpop2,fitvalue); 	
    y(i)    = bestfit;    	% ���ص� y ������Ӧ��ֵ�����Ǻ���ֵ
    
    % ���Ա��������ʮ����
    x(i,:)    = decodechrom(bestindividual, x_range, chromlength);	
    
    % ������Ⱥ
    pop     = newpop2;
end


% ��ͼ
figure(1)
plot(1:iter_num, 1./fitness_ave, 'r', 1:iter_num, 1./y, 'b')
grid on
legend('ƽ����Ӧ��', '������Ӧ��')

%Ҫ�ҳ�һ���㣬���������������е�ľ���֮�����
point = [1.4,2.7,1.5,4.6,5.2,5.6,8.2,3.8,4.6,8.7;
     	 3.6,0.1,6.9,3.6,1.2,2.7,3.5,2.1,2.9,3.3];
figure(2)
scatter(point(1,:),point(2,:), 'ko')
hold on
[z ,index]=max(y);             %�������ֵ����λ��
chrom = x(index,:);
scatter(chrom(1), chrom(2), 'bo', 'filled')
for i = 1:10
    plot([point(1, i) chrom(1)], [point(2, i) chrom(2)], 'r')
end

%%������
ymax = z;
disp(['����Ⱦɫ��Ϊ', num2str(x(index,:))])
disp(['������Ӧ��Ϊ', num2str(1./ymax)])
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
function [objvalue]=calobjvalue(pop, x_range, chromlength)
% ��֪���λ��
point = [1.4,2.7,1.5,4.6,5.2,5.6,8.2,3.8,4.6,8.7;
      	 3.6,0.1,6.9,3.6,1.2,2.7,3.5,2.1,2.9,3.3];
x1 = decodechrom(pop, x_range, chromlength);   %��popÿ��ת����ʮ������
[px,py] = size(x1);
for i = 1:px
    x = x1(i,1);
    y = x1(i,2);
    xy = [x;y] * ones(1,10);
    objvalue(i,1) =  1/sum(sqrt((xy(1,:)-point(1,:)).^2+(xy(2,:)-point(2,:)).^2));
end

end


% 2.2.2 �������Ʊ���ת��Ϊʮ������(2)
% decodechrom.m�����Ĺ����ǽ�Ⱦɫ��(������Ʊ���)ת��Ϊʮ���ƣ�
% ����pop ��ʾ��Ⱥ��
% x_range ������ȡֵ��Χ
% chromlength  �����Ʊ���ĳ���
% (���ڶ���������ԣ�������������������2*chromlengthλ��ʾ��ÿ������chromlengthλ
% ���һ��������1��ʼ����һ��������chromlength+1��ʼ)
function pop2=decodechrom(pop,x_range,chromlength)          %1  10
[num_var, temp] = size(x_range);

if num_var == 1         % ���ݱ�����ȡֵ��Χ�ж��Ƿ��ж��������num_varΪ1����ʾֻ��һ������
    pop1 = pop;
    x_min = x_range(1);
    x_max = x_range(2);    
    pop2 = x_min + (x_max - x_min) * decodebinary(pop1)/(2^chromlength - 1);  %��popÿ��ת����ʮ����ֵ
elseif num_var == 2    % 2������ʱ
    [px,py] = size(pop);
    pop1 = pop;
    popx = pop1(1:px,1:chromlength);
    popy = pop1(1:px,chromlength+1:2*chromlength);
    x_min = x_range(1);
    x_max = x_range(3);
    pop2(1:px,1) = x_min + (x_max - x_min) * decodebinary(popx)/(2^chromlength - 1);
    pop2(1:px,2) = x_min + (x_max - x_min) * decodebinary(popy)/(2^chromlength - 1);    
end

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
function fitvalue=calfitvalue(objvalue)

[px,py]=size(objvalue);                   % Ŀ��ֵ�����и�
for i=1:px
    if objvalue(i)>0
        temp=objvalue(i);
    else
        temp=0.0;
    end
    fitvalue(i)=temp;
end
fitvalue=fitvalue';

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
function [newpop]=crossover(pop,pc)                 % pc=0.6
[px,py]=size(pop);
newpop=ones(size(pop));
for i=1:2:px-1                                   	% ����Ϊ2���ǽ����ڵ�����������н���
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
