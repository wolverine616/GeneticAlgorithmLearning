%% 运动规划课程遗传算法仿真实验
% 功能： 
% 给定已知点序列，寻找到其他点位置最近的点
% Created by Yuanlong Zhang
% Date: 2020.10.19

function GeneticAlgorithmSiumlation3

clear;
clc;                        % 清屏
close all                   % 关闭所有画图窗口
x_max   = 10;               % 自变量x 的最大值
x_min   = -10;              % 自变量x 的最小值
x_range = [x_min, x_max;x_min, x_max];
[num_var,~] = size(x_range);
popsize = 200;           	% 群体大小
chromlength = 10;          	% 字符串长度(个体长度)，需要根据问题求解的精度、变量的取值范围综合判定
pc = 0.6;               	% 交叉概率，只有在随机数小于pc时，才会产生交叉 一般取 60~100%
pm = 0.02;                	% 变异概率，一般 0.1~10%
iter_num    = 1000;         % 遗传代数

pop = initpop(popsize, chromlength, num_var);             % 随机产生二进制编码的初始群体
for i = 1:1:iter_num
    % 计算目标函数
    objvalue    = calobjvalue(pop, x_range, chromlength);      
    
    % 计算群体中每个个体的适应度
    fitvalue    = calfitvalue(objvalue);          
    
    % 选择，复制
    newpop      = selection(pop,fitvalue);    
    
    % 交叉
    newpop1     = crossover(newpop,pc);            
    
    % 变异
    newpop2     = mutation(newpop1,pm);                 
    
    % 计算目标函数
    objvalue    = calobjvalue(newpop2, x_range, chromlength); 
    
    % 计算群体中每个个体的适应度
    fitvalue    = calfitvalue(objvalue);            
    
    % 平均适应度
    fitness_ave(i) = calavefitness(fitvalue);         	
    
    % 求出群体中适应值最大的个体及其适应值
    [bestindividual,bestfit]=best(newpop2,fitvalue); 	
    y(i)    = bestfit;    	% 返回的 y 是自适应度值，而非函数值
    
    % 将自变量解码成十进制
    x(i,:)    = decodechrom(bestindividual, x_range, chromlength);	
    
    % 更新种群
    pop     = newpop2;
end


% 作图
figure(1)
plot(1:iter_num, 1./fitness_ave, 'r', 1:iter_num, 1./y, 'b')
grid on
legend('平均适应度', '最优适应度')

%要找出一个点，这个点距离其他所有点的距离之和最短
point = [1.4,2.7,1.5,4.6,5.2,5.6,8.2,3.8,4.6,8.7;
     	 3.6,0.1,6.9,3.6,1.2,2.7,3.5,2.1,2.9,3.3];
figure(2)
scatter(point(1,:),point(2,:), 'ko')
hold on
[z ,index]=max(y);             %计算最大值及其位置
chrom = x(index,:);
scatter(chrom(1), chrom(2), 'bo', 'filled')
for i = 1:10
    plot([point(1, i) chrom(1)], [point(2, i) chrom(2)], 'r')
end

%%输出结果
ymax = z;
disp(['最优染色体为', num2str(x(index,:))])
disp(['最优适应度为', num2str(1./ymax)])
end


%% 2.1初始化(编码)
% initpop.m函数的功能是实现群体的初始化，popsize表示群体的大小，chromlength表示染色体的长度(二值数的长度)，
% rand随机产生每个单元为 {0,1} 行数为popsize，列数为chromlength的矩阵，
% round对矩阵的每个单元进行圆整四舍五入。
% 长度大小取决于变量的二进制编码的长度*变量个数。
% 初始化
function pop=initpop(popsize,chromlength, num_var) 
if num_var == 1                 % 只有一个变量的情况
    pop = round(rand(popsize,chromlength)); 
else                            % 多个变量的情况
    pop = round(rand(popsize,chromlength*num_var));
end

end

%% 2.2 计算目标函数值
% 2.2.1 calobjvalue.m函数的功能是实现目标函数的计算
% 实现目标函数的计算，将 二值域 中的数转化为 变量域的数
% 参数pop 表示种群数
% x_range 变量的取值范围
% chromlength  二进制编码的长度
function [objvalue]=calobjvalue(pop, x_range, chromlength)
% 已知点的位置
point = [1.4,2.7,1.5,4.6,5.2,5.6,8.2,3.8,4.6,8.7;
      	 3.6,0.1,6.9,3.6,1.2,2.7,3.5,2.1,2.9,3.3];
x1 = decodechrom(pop, x_range, chromlength);   %将pop每行转化成十进制数
[px,py] = size(x1);
for i = 1:px
    x = x1(i,1);
    y = x1(i,2);
    xy = [x;y] * ones(1,10);
    objvalue(i,1) =  1/sum(sqrt((xy(1,:)-point(1,:)).^2+(xy(2,:)-point(2,:)).^2));
end

end


% 2.2.2 将二进制编码转化为十进制数(2)
% decodechrom.m函数的功能是将染色体(或二进制编码)转换为十进制，
% 参数pop 表示种群数
% x_range 变量的取值范围
% chromlength  二进制编码的长度
% (对于多个变量而言，如有两个变量，采用2*chromlength位表示，每个变量chromlength位
% 则第一个变量从1开始，另一个变量从chromlength+1开始)
function pop2=decodechrom(pop,x_range,chromlength)          %1  10
[num_var, temp] = size(x_range);

if num_var == 1         % 根据变量的取值范围判定是否有多个变量；num_var为1，表示只有一个变量
    pop1 = pop;
    x_min = x_range(1);
    x_max = x_range(2);    
    pop2 = x_min + (x_max - x_min) * decodebinary(pop1)/(2^chromlength - 1);  %将pop每行转换成十进制值
elseif num_var == 2    % 2个变量时
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

% 将二进制数转化为十进制数
function pop2=decodebinary(pop)
[px,py]=size(pop);                  % 求pop行和列数
for i=1:py
    pop1(:,i)=2.^(py-i).*pop(:,i);  % 产生 [2^n 2^(n-1) ... 1] 的行向量，然后求和
end
pop2=sum(pop1,2);                   % 求pop1的每行之和，得到对应的十进制数

end

%% 2.3 计算个体的适应值
% 遗传算法子程序
function fitvalue=calfitvalue(objvalue)

[px,py]=size(objvalue);                   % 目标值有正有负
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

%% 求平均适应度
function fitness_ave = calavefitness(fitness)
[N ,~] = size(fitness);
fitness_ave = sum(fitness)/N;

end


%% 2.4 选择
% 选择复制
% 选择或复制操作是决定哪些个体可以进入下一代。程序中采用赌轮盘选择法选择，这种方法较易实现。
% 根据方程 pi=fi/∑fi=fi/fsum ，选择步骤：
% 1） 在第 t 代，由（1）式计算 fsum 和 pi 
% 2） 产生 {0,1} 的随机数 rand( .)，求 s=rand( .)*fsum
% 3） 求 所有fi≥s 中最小的 k ，则第 k 个个体被选中
% 4） 进行 N 次2）、3）操作，得到 N 个个体，成为第 t=t+1 代种群
% 遗传算法子程序
function [newpop]=selection(pop, fitvalue) 
totalfit = sum(fitvalue);                       % 求适应值之和
fitvalue = fitvalue/totalfit;                   % 单个个体被选择的概率
fitvalue = cumsum(fitvalue);                    % 如 fitvalue=[1 2 3 4]，则cumsum(fitvalue)=[1 3 6 10],要累加，轮盘赌法，依次看是否在转得的区域内 
[row_p,col_p] = size(pop);                      % row_p*col_p
ms = sort(rand(row_p,1));                       % 从小到大排列
fitin = 1;
newin = 1;
while newin <= row_p                           	% 选出col_p个新个体，有重复情况，和上面介绍的方法不太一样
    if(ms(newin))<fitvalue(fitin)
        newpop(newin,:) = pop(fitin,:);
        newin = newin + 1;
    else
        fitin = fitin + 1;
    end
end

end


%% 2.5 交叉
% 交叉(crossover)，群体中的每个个体之间都以一定的概率 pc 交叉，即两个个体从各自字符串的某一位置
% （一般是随机确定）开始互相交换，这类似生物进化过程中的基因分裂与重组。例如，假设2个父代个体x1，x2为：
% x1=0100110
% x2=1010001
% 从每个个体的第3位开始交叉，交又后得到2个新的子代个体y1，y2分别为：
% y1＝0100001
% y2＝1010110
% 这样2个子代个体就分别具有了2个父代个体的某些特征。利用交又我们有可能由父代个体在子代组合成具有更高适合度的个体。
% 事实上交叉是遗传算法区别于其它传统优化方法的主要特点之一。
% 遗传算法子程序
function [newpop]=crossover(pop,pc)                 % pc=0.6
[px,py]=size(pop);
newpop=ones(size(pop));
for i=1:2:px-1                                   	% 步长为2，是将相邻的两个个体进行交叉
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

%% 2.6 变异
% 变异(mutation)，基因的突变普遍存在于生物的进化过程中。变异是指父代中的每个个体的每一位都以概率 pm 翻转，
% 即由“1”变为“0”，或由“0”变为“1”。
% 遗传算法的变异特性可以使求解过程随机地搜索到解可能存在的整个空间，因此可以 在一定程度上 求得全局最优解。
% 遗传算法子程序
function [newpop]=mutation(pop,pm)
[px,py]=size(pop);
newpop=ones(size(pop));
for i=1:px
    if(rand<pm)
        mpoint=round(rand*py);              % 产生的变异点在1-10之间
        if mpoint<=0
            mpoint=1;                   	% 变异位置
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


%% 2.7 求出群体中最大得适应值及其个体
% 遗传算法子程序
% 求出第 t 代群体中适应值最大的值
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
