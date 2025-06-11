%% write your code here
clc; clear; close all;
clear strategy_Thermostat_L % 清空策略函数内部的 persistent 类型的变量

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
TAsset = TAsset(TAsset.Time >= startDate_filter & TAsset.Time < endDate_filter, :);

if isempty(TAsset)
    error('在指定日期范围 (2018-01-01 to 2025-01-01) 内没有数据，请检查CSV文件内容和日期范围。');
end

%% P2, Model 
% 2.1 分割投资的时间段：按日度分割
all_trading_days = unique(TAsset.Time);
nPeriod = length(all_trading_days);
PeriodBegins = all_trading_days;

history_window_size = 60; 

% 2.2 滚动进行投资决策
fprintf('Begin daily backtest:\n');

cellT_Portfolio = cell(nPeriod, 1);
cellT_Perform = cell(nPeriod, 1);

for iP = (history_window_size + 1):nPeriod
    % (1) 准备历史数据 和 当前数据
    current_day = PeriodBegins(iP);
    
    history_start_day = PeriodBegins(iP - history_window_size);
    history_end_day   = PeriodBegins(iP - 1);
    idx_History = (TAsset.Time >= history_start_day) & (TAsset.Time <= history_end_day);
    T_History = TAsset(idx_History, :);
    
    idx_Current = TAsset.Time == current_day;
    T_Current = TAsset(idx_Current, :);
    
    if isempty(T_Current)
        continue;
    end
    
    R_Backup = T_Current.R; 
    T_Current.R(:) = NaN;
    T_Current.P(:) = NaN;

    % (2) 构建策略函数
    TPortfolio = strategy_Thermostat_L(T_History, T_Current);
    
    % (3) 度量策略表现
    T_Current.R = R_Backup;
    T_Perform = strategy_perform(TPortfolio, T_Current);
    
    % 记录，策略的投资权重，及策略表现数据
    TPortfolio.Time(:) = current_day;
    cellT_Portfolio{iP} = TPortfolio;
    cellT_Perform{iP} = T_Perform;
end
fprintf('\nEnd at: %g, %s;\n', iP, datestr(PeriodBegins(iP)));

% 2.3 数据整理：完成循环之后，汇总策略，及其表现
T_PerformAll = vertcat(cellT_Perform{:});
T_PerformAll.RCum = cumprod(T_PerformAll.R + 1, 1);

%% P3, 投资结果展示：策略及其表现

% 3.1 准备用于对比的上证指数数据
T_Benchmark = T_SH_Index(ismember(T_SH_Index.Time, T_PerformAll.Time), :);
plot_benchmark = false;

if ~isempty(T_Benchmark)
    T_Benchmark.NormalizedClose = T_Benchmark.Close / T_Benchmark.Close(1);
    plot_benchmark = true;
else
    warning('上证指数数据与策略回测日期无法对齐，将不绘制对比曲线。');
end

% 3.2 投资策略表现绘图
hf = figure('Name', 'Thermostat策略与上证指数累计净值对比');
plot(T_PerformAll.Time, T_PerformAll.RCum, 'b-', 'LineWidth', 2);
hold on;

if plot_benchmark
    plot(T_Benchmark.Time, T_Benchmark.NormalizedClose, 'r-', 'LineWidth', 2);
    legend('Thermostat策略累计净值', '上证指数', 'Location', 'northwest');
else
    legend('Thermostat策略累计净值', 'Location', 'northwest');
end

hold off;
title('Thermostat策略与上证指数累计净值对比');
xlabel('时间');
ylabel('累计净值 (归一化)');
grid on;
datetick('x', 'yyyy-mm');

% 创建用于存放图片的文件夹 (如果不存在)
s_Path = '91_程序运行结果输出';
if ~exist(s_Path, 'dir')
    mkdir(s_Path);
end
% 定义图片文件路径及格式
s_FigPath = fullfile(s_Path, '投资策略累计回报率.emf');
% 保存图片
saveas(hf, s_FigPath);
fprintf('\n策略表现图已保存至: %s\n', s_FigPath);


%% P4, 绩效指标
FinalReturn = T_PerformAll.RCum(end);
AnnualReturn = FinalReturn^(252/height(T_PerformAll)) - 1;
MaxDrawdown = maxdrawdown(T_PerformAll.RCum);
SharpeRatio = mean(T_PerformAll.R, 'omitnan') / std(T_PerformAll.R, 'omitnan') * sqrt(252);

% 保留您要求的性能评估报告
fprintf('\n=============== Thermostat策略评估 (日度) ===============\n');
fprintf('回测区间: %s to %s\n', datestr(T_PerformAll.Time(1),'yyyy-mm-dd'), datestr(T_PerformAll.Time(end),'yyyy-mm-dd'));
fprintf('最终累计收益率: %.2f%%\n', (FinalReturn - 1) * 100);
fprintf('年化收益率: %.2f%%\n', AnnualReturn * 100);
fprintf('最大回撤: %.2f%%\n', MaxDrawdown * 100);
fprintf('夏普比率 (年化): %.2f\n', SharpeRatio);
fprintf('====================================================\n');