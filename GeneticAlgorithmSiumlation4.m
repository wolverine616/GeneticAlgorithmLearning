%% 运动规划课程遗传算法仿真实验
% 功能： 
% 给定已知点序列，寻找到其他点位置最近的点
% Created by Yuanlong Zhang
% Date: 2020.10.19

function GeneticAlgorithmSiumlation4

clear;
clc;                        % 清屏
close all                   % 关闭所有画图窗
x_max   = 20;               % 坐标的最大值，默认设置横坐标 x 和纵坐标 y 的刻度相同
x_min   = 0;                % 坐标的最小值
N_point = 10;            	% 搜索的路径点个数，可根据需要自行设定
x_range = repmat([x_min, x_max],[N_point*2,1]);    % N_point个点，2*N_point个变量， repmat（A，[m,n]),产生m*n个A
[num_var,~] = size(x_range);
popsize = 200;           	% 群体大小
chromlength = 10;          	% 字符串长度(个体长度)，需要根据问题求解的精度、变量的取值范围综合判定
pc = 0.6;               	% 交叉概率，只有在随机数小于pc时，才会产生交叉 一般取 60~100%
pm = 0.02;                	% 变异概率，一般 0.1~10%
iter_num    = 200;        	% 遗传求解的迭代数

% 首先，利用循环产生一个可行的初始路径
isTerminal  = 0;            % 循环结束标志位，其值为0或1
iter_max = 100;          	% 初始最大迭代次数，默认设置 100；
                            % 超过仍不能找到合适的初始解，表示当前搜索路径问题无可行解，可增大或减小搜索路径
Counter_i = 1;
while ~isTerminal           % ~ 取反符号；isTerminal=0，继续循环；1，结束循环
    % 编码，随机产生二进制编码的初始群体
    pop = initpop(popsize, chromlength, num_var);    
    
    % 计算目标值
    [objvalue, x, y] = calobjvalue(pop, x_range, chromlength);
    
    % 计算群体中每个个体的适应度
    fitvalue = calfitvalue(x, y, objvalue);
    
    if ~isempty(find(fitvalue)>0) || Counter_i >= iter_max     % 若不全为0，则说明找到了一个合理的初值，isTerminal =1，结束循环
         isTerminal =1;                                        % 同时，当迭代次数Counter_i 超过设定值时也必须结束循环       
    else
        Counter_i = Counter_i + 1;
    end    
end

if Counter_i >= iter_max
    error('当前无可行解！请检查问题参数设置！')
end

for i = 1:1:iter_num
    % 选择，复制
    newpop      = selection(pop,fitvalue);     
    
    % 交叉
    newpop1     = crossover(newpop,pc);                 
    
    % 变异
    newpop2     = mutation(newpop1,pm);               
    
    % 计算目标函数

    [objvalue, x, y] = calobjvalue(newpop2, x_range, chromlength);
    
    % 计算群体中每个个体的适应度
    fitvalue = calfitvalue(x, y, objvalue);
    
    % 平均适应度
    fitness_ave(i) = calavefitness(fitvalue);         	
    
    % 求出群体中适应值最大的个体及其适应值
    [bestindividual,bestfit]=best(newpop2,fitvalue); 
    
    % 返回的 best_value 是自适应度值，而非函数值
    best_value(i)    = bestfit;        
    
    % 解码，将自变量解码成十进制数
    [x_best(i,:), y_best(i,:)]    = decodechrom(bestindividual, x_range, chromlength);	
    
    % 更新种群
    pop     = newpop2;
end


% 作图
figure(1)
plot(1:iter_num, fitness_ave, 'r', 1:iter_num, best_value, 'b')
grid on
legend('平均适应度', '最优适应度')

% 构建问题模型画布
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

[z ,index]=max(best_value);                     % 计算最大值及其位置
x_chrom_ = x_best(index,:);                     % 最短路径 x 坐标
y_chrom_ = y_best(index,:);                     % 最短路径 y 坐标
scatter(x_chrom_, y_chrom_, 'bo', 'filled');    % 最短路径的散点图
plot(x_chrom_,y_chrom_, 'r','Linewidth',2)      % 连线，得到最短路径

% 输出结果
ymax = z;
disp(['最优染色体为', num2str([x_chrom_,y_chrom_])])
disp(['最短路径长度为', num2str(-(ymax-200))])
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
%这个函数执行完之后[objvalue , x, y],objvalue
%是一个200行，一列的矩阵，保存着每一行的x，y计算后的值,x,y的每一行是十个点


% 2.2.2 将二进制编码转化为十进制数(2)
% decodechrom.m函数的功能是将染色体(或二进制编码)转换为十进制，
% 参数pop 表示种群数
% x_range 变量的取值范围
% chromlength  二进制编码的长度
% (对于多个变量而言，如有两个变量，采用2*chromlength位表示，每个变量chromlength位
% 则第一个变量从1开始，另一个变量从chromlength+1开始)
function [x_sort, y_sort]=decodechrom(pop,x_range,chromlength)          
[num_var, temp] = size(x_range);    % 
[px,~] = size(pop);
if num_var == 1         % 根据变量的取值范围判定是否有多个变量；num_var为1，表示只有一个变量
    pop1 = pop;
    x_min = x_range(1);
    x_max = x_range(2);    
    pop2 = x_min + (x_max - x_min) * decodebinary(pop1)/(2^chromlength - 1);  %将pop每行转换成十进制值
else                    % 多个变量时
    pop1 = pop;
    x_min = x_range(1);
    x_max = x_range(21);
    for i = 1:num_var
        popvar = pop1(1:px,(i-1)*chromlength+1:i*chromlength);
        pop2(1:px,i) = x_min + (x_max - x_min) * decodebinary(popvar)/(2^chromlength - 1);
    end
end

% 将 pop2 中的数分别赋给横坐标 x 和纵坐标 y
x = pop2(:,1:num_var/2);
y = pop2(:,num_var/2+1:end);

% 要求初始从（0，0），终端为（20，20），因此修改初始和终端值
%此时x,y分别是一个200行，10列的矩阵。
x(:,1) = 0;
x(:,num_var/2) = 20;
y(:,1) = 0;
y(:,num_var/2) = 20;
%此时x,y矩阵的第一行和最后一行分别是0和20，十个数字一行表明十个点,正好是一条路

% 重新按大小进行排序
%行按大小排列
x_sort = sort(x,2);     % 调用 sort 函数，按行升序排列
y_sort = sort(y,2);

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
% 计算个体的适应值
function fitvalue = calfitvalue(x, y, objvalue)

%设置环境模型
%x3、y3、x4、y4 为所有障碍物边界的相邻顶点组成的线段。
x3=[2,4,4,2,2,4,4,0,0,2,2,6,6,2,0,2,2,0,6,10,10,6,8,10,10,8,12,14,14,12,10,12,12,10,12,14,14,12,14,18,18,16,16,14,16,20,20,16,18,20,20,18];
y3=[0,0,2,2,4,4,8,8,6,6,12,12,14,14,18,18,20,20,2,2,6,6,16,16,18,18,6,6,8,8,12,12,14,14,18,18,20,20,12,12,14,14,16,16,2,2,4,4,8,8,10,10];
x4=[4,4,2,2,4,4,0,0,2,2,6,6,2,2,2,2,0,0,10,10,6,6,10,10,8,8,14,14,12,12,12,12,10,10,14,14,12,12,18,18,16,16,14,14,20,20,16,16,20,20,18,18,];
y4=[0,2,2,0,4,8,8,6,6,4,12,14,14,12,18,20,20,18,2,6,6,2,16,18,18,16,6,8,8,6,12,14,14,12,18,20,20,18,12,14,14,16,16,12,2,4,4,2,8,10,10,8];

[row_x,col_x]=size(x);          % 取 x 的行和列
fitvalue = zeros(row_x,1);      % 初始化列向量 fitvalue
%row_x=200 , col_x = 10;
for j=1:1:row_x
    % （x1,y1)存放第j行的第1到倒数第二个点
    x1 = x(j,1:col_x-1);
    y1 = y(j,1:col_x-1);
    
    
    % （x2,y2)存放第j行的第2到倒数第一个点
    x2 = x(j,2:col_x);
    y2 = y(j,2:col_x);
    
    %x1,y1一行一条路十个点的前九个点的矩阵
    %x2,y2一行一条路十个点的后九个点的矩阵
    % 判断种群相邻点连成的线段是否与所有障碍物边界相交
    ch = check(x1,y1,x2,y2,x3,y3,x4,y4);   
    
    % 如果矩阵 [ch] 所有元素的和不等于0,则说明存在某一条线段经过障碍物边界。
    if sum(sum(ch))>0           
        temp = 0;
    else
        temp = objvalue(j);
    end
    fitvalue(j)=temp;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%函数名称：边界相角检测函数 check.m
%%入口参数：两条线段的顶点。例如 判断线段PQ 与 XY 是否相交，其中 P（x1,y1）、Q(x2,y2)、X(x3,y3)、Y(x4,y4)
                            %%则函数参数设置为：chack(x1,y1,x2,y2,x3,y3,x4,y4)
%%出口参数：　是否相交 0：不相交  1：相交
%%说明：
    %%该函数可以检测两条线断是否相交，当入口参数是向量的时候，可以同时检测多条线段是否相交。
    %%该函数入口参数可以为向量，向量可以是列向量也可以是行向量，但必须保证所有入口参数的维度相同。
    %%检测方式采用相互跨立实验实现。当两条线段相互跨立时，则两线段在二维平面中相交。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pop]=check(x1,y1,x2,y2,x3,y3,x4,y4)
%x1,y1一行一条路十个点的前九个点的矩阵
    %x2,y2一行一条路十个点的后九个点的矩阵
[px,py]=size(x1);
p=max(px,py);       % 实现行向量列向量都可以进行检测
[px,py]=size(x3);
q=max(px,py);
%此时px=1,py=9;p=9;
%px=1,py=52;q=52
%采用相互跨立实验
for j=1:1:q         % 检测线段相交的次数
    for i=1:1:p     % 检测所有种群相邻点连成的线段是否与当前障碍物边界相交
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
        
        
        if m*n<=0 && mm*nn<=0     % 相交
            pop(j,i)=1;           % 将结果存到矩阵[pop]中。
        else                      % 显然矩阵[pop]是一个q行p列的矩阵。
            pop(j,i)=0;
        end
    end
end
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
function [newpop]=crossover(pop,pc)                         % pc=0.6
[px,py]=size(pop);
newpop=ones(size(pop));
for i=1:2:px-1                                              % 步长为2，是将相邻的两个个体进行交叉
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

