% clear all
% mex cec17_func.cpp -DWINDOWS

fhd=str2func('cec17_func');
%func_num =17
Popsize=100; 
%Z=123;
dim=30;
lb = -100.*ones(1,dim);    % 取值下限向量
ub = 100.*ones(1,dim);    % 取值上限向量
%ub=100;
%lb=-100;
Multi_Num=1;
Max_iteration=10000;
vmax = 0.2 * (ub - lb); % vmax为搜索区间的20%(PSO)
vmin = -vmax; % vmin为-vmax

data_bdbo=[];
bdbo_Wilcoxon=[];

for i =20:26
    func_num=i;
    bestfit_BDBO_Z=[];
    CEC=i
    MCBDBO=[];

    for j=123:133
        Z=j;   
        
        %% BDBO收敛曲线
        [Multi_Convergence_curve_BDBO, Multi_time, bestXfit_BDBO, GMM_EM] = bdbo(Z,Popsize, dim, ub, lb, @fobj, Max_iteration, Multi_Num,func_num);
        bestfit_BDBO_Z=[bestfit_BDBO_Z;bestXfit_BDBO];
        MCBDBO=[MCBDBO,Multi_Convergence_curve_BDBO];
        %% BDBO迭代数据
        %data_bdbo=[data_bdbo;[func_num,min(bestfit_BDBO_Z),mean(bestfit_BDBO_Z),var(bestfit_BDBO_Z)]];      

    end
    
    figure;
    
    % 绘制DBO的收敛曲线
    %plot(mean(MCDBO,2), 'LineWidth', 2);%DBO
    %hold on; % 使得下一次plot在同一张图上

    %hold on; % 使得下一次plot在同一张图上
    plot(mean(MCBDBO,2), 'LineWidth', 2);%BDBO

    % 添加图例、标题和轴标签
    legend('BDBO')%
    title('Convergence Curves');
    title(['Convergence Curves for Function ',',func_num  = ' num2str(func_num)]);
    xlabel('Iteration');
    ylabel('Value');
        
    % 添加网格线
    %grid on;

    data_bdbo=[data_bdbo;[CEC,min(bestX_Z),mean(bestX_Z),std(bestX_Z)]];
    bdbo_Wilcoxon=[bdbo_Wilcoxon;[CEC,bestX_Z(:)']];
    
    
end
format longG;

data_bdbo

bdbo_Wilcoxon
p_Wilcoxon=[];
h_Wilcoxon=[];
for i=1:size(gdbo_Wilcoxon,1)
    [pbdbo, hbdbo] = ranksum(bdbo_Wilcoxon(i,2:end), dbo_Wilcoxon(i,2:end)); 
    p_Wilcoxon=[p_Wilcoxon;[gdbo_Wilcoxon(i,1),pbdbo,prdbo,pldbo,pgdbo]];
    h_Wilcoxon=[h_Wilcoxon;[gdbo_Wilcoxon(i,1),hbdbo,hrdbo,hldbo,hgdbo]];
end
p_Wilcoxon
h_Wilcoxon



t=1;
% 创建一个新的图形窗口
%figure; % 创建一个新的图形窗口
%plot(Multi_Convergence_curve); % 绘制收敛曲线
%title('Convergence Curve'); % 给图形添加标题
%xlabel('Iteration'); % 给x轴添加标签
%ylabel('Value'); % 给y轴添加标签
%grid on; % 添加网格线
