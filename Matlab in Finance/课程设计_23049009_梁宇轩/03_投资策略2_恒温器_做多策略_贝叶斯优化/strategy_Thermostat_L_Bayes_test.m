%% write your code here
clc; clear; close all;
clear strategy_Thermostat_L objective_Thermostat
global TAsset_global; % 定义全局变量，用于向目标函数传递数据

%% 目录路径处理
p1 = fileparts(which(mfilename));
addpath(genpath(p1));    % 将当前目录和子目录，添加到搜索路径
cd(p1)

%% Part 1. inputs
% 读取干净的原始数据
dataPath = fullfile(p1, '01_数据及其函数');
[TAsset, T_SH_Index] = Data_Import_Thermostat(dataPath); 

% 初始日期筛选
startDate_filter = datetime('2018-01-01');
endDate_filter = datetime('2025-01-01');
TAsset_global = TAsset(TAsset.Time >= startDate_filter & TAsset.Time < endDate_filter, :);

if isempty(TAsset_global)
    error('在指定日期范围 (2018-01-01 to 2025-01-01) 内没有数据。');
end

%% P2, Model 
% 此部分为贝叶斯优化过程，用于寻找最优模型参数
fprintf('====== 开始贝叶斯优化 ======\n');

% 2.1 定义待优化的参数及其范围
optimVars = [
    optimizableVariable('swingTrendSwitch', [10, 50], 'Type', 'integer')
    optimizableVariable('swingPrcnt1', [0.2, 1.0], 'Type', 'real')
    optimizableVariable('bollingerLengths', [20, 100], 'Type', 'integer')
    optimizableVariable('numStdDevs', [1.5, 3.5], 'Type', 'real')
    optimizableVariable('atrLength', [5, 30], 'Type', 'integer')
];

% 2.2 执行贝叶斯优化
% 优化时，动态显示'最小目标函数值'的收敛过程，保持界面简洁
BayesObject = bayesopt(@objective_Thermostat, optimVars, ...
    'MaxObjectiveEvaluations', 30, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'PlotFcn', @plotMinObjective, ...
    'Verbose', 1);

% 2.3 提取最优参数
best_params = BayesObject.XAtMinObjective;
fprintf('====== 贝叶斯优化完成 ======\n');
disp('找到的最优参数:');
disp(best_params);


%% P3, 回测与可视化
fprintf('\n====== 使用最优参数进行回测 ======\n');

% 3.1 使用最优参数进行回测
[T_PerformAll, T_Turnover] = backtest_Thermostat_L(TAsset_global, best_params);

% 3.2 准备基准指数数据
plot_benchmark = false;
if exist('T_SH_Index','var') && ~isempty(T_SH_Index)
    T_Benchmark = T_SH_Index(T_SH_Index.Time >= T_PerformAll.Time(1) & T_SH_Index.Time <= T_PerformAll.Time(end),:);
    T_Benchmark.NormalizedClose = T_Benchmark.Close / T_Benchmark.Close(1);
    plot_benchmark = true;
end

% 3.3 投资策略表现绘图
hf_perf = figure('Name', '最优参数下的策略表现');
plot(T_PerformAll.Time, T_PerformAll.RCum, 'b-', 'LineWidth', 2);
hold on;
if plot_benchmark
    plot(T_Benchmark.Time, T_Benchmark.NormalizedClose, 'r-', 'LineWidth', 2);
    legend('优化后策略净值', '上证指数', 'Location', 'northwest');
else
    legend('优化后策略净值', 'Location', 'northwest');
end
hold off;
title('最优参数下的Thermostat策略表现');
xlabel('时间'); ylabel('累计净值 (归一化)');
grid on; datetick('x', 'yyyy-mm');

% 3.4 准备保存路径
s_Path = '91_程序运行结果输出';
if ~exist(s_Path, 'dir'), mkdir(s_Path); end

% 3.5 导出回测表现图
s_FigPath_Perf = fullfile(s_Path, '投资策略累计回报率_优化后.emf');
saveas(hf_perf, s_FigPath_Perf);
fprintf('\n策略表现图已保存至: %s\n', s_FigPath_Perf);


% --- 代码修改点: 新增的参数收敛过程可视化 ---
% 3.6 绘制并保存贝叶斯优化参数收敛过程图
hf_convergence = figure('Name', '贝叶斯优化参数收敛过程', 'Position', [100, 100, 800, 1000]);
paramTrace = BayesObject.XTrace; % 提取所有迭代的参数值
iterations = 1:height(paramTrace);
paramNames = paramTrace.Properties.VariableNames;
bestParamsTable = BayesObject.XAtMinObjective;

for i = 1:length(paramNames)
    subplot(length(paramNames), 1, i);
    plot(iterations, paramTrace.(paramNames{i}), '.-', 'MarkerSize', 12);
    hold on;
    % 用红色虚线标出最终找到的最优值
    yline(bestParamsTable.(paramNames{i}), 'r--', 'LineWidth', 1.5, {'最终最优值'});
    hold off;
    title(paramNames{i}, 'Interpreter', 'none'); % 'Interpreter', 'none' 防止变量名中的下划线被识别为下标
    ylabel('参数值');
    grid on;
    if i < length(paramNames)
        set(gca, 'XTickLabel', []); % 隐藏非底部子图的x轴标签
    end
end
xlabel('迭代次数');
sgtitle('贝叶斯优化参数收敛过程', 'FontSize', 14, 'FontWeight', 'bold'); % 为整个图窗添加总标题

% 3.7 导出参数收敛图
s_FigPath_Conv = fullfile(s_Path, '贝叶斯参数收敛过程图.emf');
saveas(hf_convergence, s_FigPath_Conv);
fprintf('贝叶斯参数收敛过程图已保存至: %s\n', s_FigPath_Conv);
% --- 修改结束 ---


%% P4, 绩效指标
FinalReturn = T_PerformAll.RCum(end);
AnnualReturn = FinalReturn^(252/height(T_PerformAll)) - 1;
AnnualVolatility = std(T_PerformAll.R) * sqrt(252);
SharpeRatio = AnnualReturn / AnnualVolatility;
MaxDrawdown = maxdrawdown(T_PerformAll.RCum);
TurnoverRate = mean(T_Turnover.Turnover);

fprintf('\n====== 策略绩效指标 ======\n');
fprintf('最终累计净值: %.4f\n', FinalReturn);
fprintf('年化收益率: %.2f%%\n', AnnualReturn * 100);
fprintf('年化波动率: %.2f%%\n', AnnualVolatility * 100);
fprintf('夏普比率: %.4f\n', SharpeRatio);
fprintf('最大回撤: %.2f%%\n', MaxDrawdown * 100);
fprintf('平均换手率: %.4f\n', TurnoverRate);
fprintf('========================\n');