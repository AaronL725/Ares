%% =======================================================================
% 量化策略参数对比可视化分析
% =======================================================================
clc; clear; close all;

%% 路径管理
p1 = fileparts(which(mfilename));
addpath(genpath(p1));
cd(p1);

% 创建输出文件夹
s_Path = '91_程序运行结果输出';
if ~exist(s_Path, 'dir'), mkdir(s_Path); end
fprintf('所有图表和结果将保存在文件夹: %s\n', s_Path);

%% =======================================================================
% Part 1: 定义参数与加载数据
% =======================================================================
% --- 策略参数组 ---
params1.name = '贝叶斯优化参数';
params1.swingTrendSwitch = 19;
params1.swingPrcnt1 = 0.76205;
params1.bollingerLengths = 31;
params1.numStdDevs = 2.6911;
params1.cmi_period = 30;
params1.atrLength = 20;

params2.name = '普通参数';
params2.swingTrendSwitch = 20;
params2.swingPrcnt1 = 10.50;
params2.bollingerLengths = 50;
params2.numStdDevs = 2;
params2.cmi_period = 30;
params2.atrLength = 10;

% --- 数据加载 ---
dataPath = fullfile(p1, '01_数据及其函数');
[TAsset, T_SH_Index] = Data_Import_Thermostat(dataPath);

% --- 回测区间 ---
startDate_filter = datetime('2018-01-01');
endDate_filter = datetime('2025-01-01');
TAsset = TAsset(TAsset.Time >= startDate_filter & TAsset.Time < endDate_filter, :);
T_SH_Index = T_SH_Index(T_SH_Index.Time >= startDate_filter & T_SH_Index.Time < endDate_filter, :);

if isempty(TAsset)
    error('在指定日期范围 (2018-01-01 to 2025-01-01) 内没有数据，请检查CSV文件内容和日期范围。');
end

%% =======================================================================
% Part 2: 执行回测
% =======================================================================
fprintf('\n开始执行回测...\n');
fprintf('正在运行策略: %s\n', params1.name);
[T_Perform1, TradeLog1] = run_backtest(TAsset, @strategy_Thermostat_L, params1);

fprintf('正在运行策略: %s\n', params2.name);
[T_Perform2, TradeLog2] = run_backtest(TAsset, @strategy_Thermostat_L, params2);
fprintf('回测完成。\n');

%% =======================================================================
% Part 3: 计算并输出关键性能指标 (KPI)
% =======================================================================
metrics1 = calculate_metrics(T_Perform1, ['恒温器策略 (' params1.name ')']);
metrics2 = calculate_metrics(T_Perform2, ['恒温器策略 (' params2.name ')']);
metrics_benchmark = calculate_metrics(T_SH_Index, '上证指数');

fprintf('\n--- 策略性能指标 (回测区间: %s to %s) ---\n', ...
    string(startDate_filter), string(TAsset.Time(end)));
disp(struct2table([metrics1, metrics2, metrics_benchmark]));


%% =======================================================================
% Part 4: 生成对比可视化图表
% =======================================================================
fprintf('\n开始生成可视化图表...\n');
plot_cumulative_returns(T_Perform1, T_Perform2, T_SH_Index, s_Path);
plot_radar_performance(metrics1, metrics2, s_Path);
plot_cumulative_alpha(T_Perform1, T_Perform2, s_Path);
plot_trade_density(TradeLog1, TAsset, s_Path);
plot_trade_density(TradeLog2, TAsset, s_Path);
plot_rolling_efficiency(TradeLog1, TradeLog2, s_Path, 30);
fprintf('所有图表已生成并保存。\n');

%% =======================================================================
% Part 5: 本地辅助函数 (Local Functions)
% =======================================================================

function [T_PerformAll, T_TradeLog] = run_backtest(TAsset, strategy_func, params)
    all_trading_days = unique(TAsset.Time);
    nPeriod = length(all_trading_days);
    history_window_size = 60;
    
    TPortfolio_All = cell(nPeriod, 1);
    TCurrent_All = cell(nPeriod, 1);
    
    for t = 1:nPeriod
        current_date = all_trading_days(t);
        history_start_date = all_trading_days(max(1, t - history_window_size));
        THistory = TAsset(TAsset.Time >= history_start_date & TAsset.Time < current_date, :);
        TCurrent = TAsset(TAsset.Time == current_date, :);
        
        if isempty(TCurrent) || isempty(THistory)
            continue;
        end
        
        TPortfolio = strategy_func(THistory, TCurrent, params);
        
        TPortfolio_All{t} = TPortfolio;
        TCurrent_All{t} = TCurrent;
    end
    
    TPortfolio_All = vertcat(TPortfolio_All{:});
    TCurrent_All = vertcat(TCurrent_All{:});
    TPortfolio_All.Properties.Description = params.name;

    T_PerformAll = strategy_perform(TPortfolio_All, TCurrent_All);
    T_TradeLog = generate_trade_log(TPortfolio_All, TAsset);
end

function metrics = calculate_metrics(T_Perf, name)
    if isempty(T_Perf) || height(T_Perf) < 2
        metrics = struct('Name', string(name), 'CumulativeReturn', 0, 'AnnualizedReturn', 0, 'MaxDrawdown', 0, 'SharpeRatio', NaN);
        return;
    end
    
    if ismember('RCum', T_Perf.Properties.VariableNames)
        returns = T_Perf.R;
        cum_returns = T_Perf.RCum;
        time_vector = T_Perf.Time;
    else
        if ismember('Datetime', T_Perf.Properties.VariableNames)
            T_Perf = renamevars(T_Perf, 'Datetime', 'Time');
        end
        time_vector = T_Perf.Time;
        T_Perf.R = [0; diff(T_Perf.Close) ./ T_Perf.Close(1:end-1)];
        T_Perf.RCum = cumprod(1 + T_Perf.R) - 1;
        returns = T_Perf.R;
        cum_returns = T_Perf.RCum;
    end
    
    metrics.Name = string(name);
    metrics.CumulativeReturn = cum_returns(end) * 100;
    
    annual_factor = 252 / mean(days(diff(time_vector)));
    metrics.AnnualizedReturn = mean(returns, 'omitnan') * annual_factor * 100;
    
    high_water_mark = cummax(cum_returns + 1);
    drawdown = (high_water_mark - (cum_returns + 1)) ./ high_water_mark;
    metrics.MaxDrawdown = max(drawdown) * 100;

    metrics.SharpeRatio = (mean(returns, 'omitnan') / std(returns, 'omitnan')) * sqrt(annual_factor);
end

function T_TradeLog = generate_trade_log(TPortfolio, TAsset)
    all_asset_ids = unique(TPortfolio.AssetID);
    strat_name = TPortfolio.Properties.Description;
    trades = [];
    
    for i = 1:length(all_asset_ids)
        asset_id = all_asset_ids{i};
        T_asset_prices = TAsset(strcmp(TAsset.AssetID, asset_id), :);
        T_asset_weights = TPortfolio(strcmp(TPortfolio.AssetID, asset_id), :);
        
        T_merged = outerjoin(T_asset_prices, T_asset_weights, 'Keys', {'Time', 'AssetID'}, 'MergeKeys', true, 'Type', 'left');
        T_merged.Weight = fillmissing(T_merged.Weight, 'previous', 'EndValues', 0);
        
        w_prev = [0; T_merged.Weight(1:end-1)];
        trade_signal = T_merged.Weight - w_prev;
        trade_indices = find(trade_signal ~= 0);
        
        direction = 0;
        entry_date = NaT;
        entry_price = NaN;
        
        for k = 1:length(trade_indices)
            idx = trade_indices(k);
            current_date = T_merged.Time(idx);
            current_price = T_merged.Close(idx);

            if direction == 0
                if T_merged.Weight(idx) == 1
                    direction = 1;
                    entry_date = current_date;
                    entry_price = current_price;
                elseif T_merged.Weight(idx) == -1
                    direction = -1;
                    entry_date = current_date;
                    entry_price = current_price;
                end
            elseif T_merged.Weight(idx) == 0
                pnl = (current_price - entry_price) * direction;
                trade_dir_str = '多头';
                if direction == -1, trade_dir_str = '空头'; end
                
                trades = [trades; {strat_name, asset_id, entry_date, current_date, entry_price, current_price, trade_dir_str, pnl}];
                
                direction = 0; entry_date = NaT; entry_price = NaN;
            end
        end
    end
    
    if isempty(trades)
        T_TradeLog = table();
        T_TradeLog.Properties.Description = strat_name;
        warning('策略“%s”没有产生任何交易。', strat_name);
        return;
    end
    
    T_TradeLog = cell2table(trades, 'VariableNames', {'Strategy', 'AssetID', 'EntryDate', 'ExitDate', 'EntryPrice', 'ExitPrice', 'Direction', 'PnL_per_Unit'});
    T_TradeLog.HoldingPeriod = days(T_TradeLog.ExitDate - T_TradeLog.EntryDate);
end

function plot_cumulative_returns(T_Perf1, T_Perf2, T_Benchmark_orig, save_path)
    T_Benchmark = T_Benchmark_orig;
    hf = figure('Name', '累计收益对比图', 'Position', [100, 100, 1000, 600]);
    hold on;

    if ismember('Datetime', T_Benchmark.Properties.VariableNames)
        T_Benchmark = renamevars(T_Benchmark, 'Datetime', 'Time');
    end
    T_Benchmark.R = [0; diff(T_Benchmark.Close) ./ T_Benchmark.Close(1:end-1)];
    T_Benchmark.RCum = cumprod(1 + T_Benchmark.R) - 1;

    plot(T_Perf1.Time, T_Perf1.RCum, 'b-', 'LineWidth', 2, 'DisplayName', T_Perf1.Properties.Description);
    plot(T_Perf2.Time, T_Perf2.RCum, 'g-', 'LineWidth', 2, 'DisplayName', T_Perf2.Properties.Description);
    plot(T_Benchmark.Time, T_Benchmark.RCum, 'k--', 'LineWidth', 1.5, 'DisplayName', '上证指数');

    annotate_curve(T_Perf1.Time, T_Perf1.RCum, 'b');
    annotate_curve(T_Perf2.Time, T_Perf2.RCum, 'g');
    annotate_curve(T_Benchmark.Time, T_Benchmark.RCum, 'k');

    hold off;
    title('策略累计收益率对比');
    xlabel('时间');
    ylabel('累计收益率');
    legend('show', 'Location', 'northwest');
    grid on;
    datetick('x', 'yyyy-mm', 'keeplimits');

    saveas(hf, fullfile(save_path, 'P1_Cumulative_Returns_Combined.emf'));
    fprintf('图1: 累计收益曲线图 (策略对比) -> 已保存\n');

    % --- 生成并保存 原始策略 vs 上证指数 ---
    hf_single1 = figure('Name', [T_Perf2.Properties.Description ' vs 上证指数'], 'Position', [100, 100, 1000, 600]);
    hold on;
    plot(T_Perf2.Time, T_Perf2.RCum, 'g-', 'LineWidth', 2, 'DisplayName', T_Perf2.Properties.Description);
    plot(T_Benchmark.Time, T_Benchmark.RCum, 'k--', 'LineWidth', 1.5, 'DisplayName', '上证指数');
    annotate_curve(T_Perf2.Time, T_Perf2.RCum, 'g');
    annotate_curve(T_Benchmark.Time, T_Benchmark.RCum, 'k');
    hold off;
    title([T_Perf2.Properties.Description ' 累计收益率 vs 上证指数']);
    xlabel('时间');
    ylabel('累计收益率');
    legend('show', 'Location', 'northwest');
    grid on;
    datetick('x', 'yyyy-mm', 'keeplimits');
    saveas(hf_single1, fullfile(save_path, 'P1a_Cumulative_Returns_Normal_vs_SH.emf'));
    fprintf('图1a: 累计收益曲线图 (%s vs 上证指数) -> 已保存\n', T_Perf2.Properties.Description);

    % --- 生成并保存 优化策略 vs 上证指数 ---
    hf_single2 = figure('Name', [T_Perf1.Properties.Description ' vs 上证指数'], 'Position', [100, 100, 1000, 600]);
    hold on;
    plot(T_Perf1.Time, T_Perf1.RCum, 'b-', 'LineWidth', 2, 'DisplayName', T_Perf1.Properties.Description);
    plot(T_Benchmark.Time, T_Benchmark.RCum, 'k--', 'LineWidth', 1.5, 'DisplayName', '上证指数');
    annotate_curve(T_Perf1.Time, T_Perf1.RCum, 'b');
    annotate_curve(T_Benchmark.Time, T_Benchmark.RCum, 'k');
    hold off;
    title([T_Perf1.Properties.Description ' 累计收益率 vs 上证指数']);
    xlabel('时间');
    ylabel('累计收益率');
    legend('show', 'Location', 'northwest');
    grid on;
    datetick('x', 'yyyy-mm', 'keeplimits');
    saveas(hf_single2, fullfile(save_path, 'P1b_Cumulative_Returns_Optimized_vs_SH.emf'));
    fprintf('图1b: 累计收益曲线图 (%s vs 上证指数) -> 已保存\n', T_Perf1.Properties.Description);

    function annotate_curve(time, cum_returns, color)
        if isempty(time), return; end
        plot(time(end), cum_returns(end), 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
        text(time(end), cum_returns(end), sprintf(' %.2f%%', cum_returns(end)*100), 'VerticalAlignment', 'bottom');

        high_water_mark = cummax(cum_returns + 1);
        drawdown = (high_water_mark - (cum_returns + 1)) ./ high_water_mark;
        if isempty(drawdown) || all(isnan(drawdown)) || max(drawdown) == 0, return; end
        if any(~isnan(drawdown))
            max_dd = max(drawdown(~isnan(drawdown)));
            idx_dd_end = find(drawdown == max_dd, 1, 'last');
            idx_dd_start = find((cum_returns + 1) == high_water_mark(idx_dd_end), 1, 'first');

            if ~isempty(idx_dd_start) && ~isempty(idx_dd_end)
                plot([time(idx_dd_start), time(idx_dd_end)], [cum_returns(idx_dd_start), cum_returns(idx_dd_end)], 'v-', 'Color', color, 'MarkerFaceColor', 'r', 'LineWidth', 1, 'MarkerSize', 6);
                text(time(idx_dd_end), cum_returns(idx_dd_end), sprintf(' 最大回撤: %.2f%%', max_dd*100), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', 8);
            end
        end
    end
end

function plot_radar_performance(metrics1, metrics2, save_path)
    labels = {'累计回报率', '年化回报率', '夏普比率', '逆最大回撤'};
    
    data1 = [metrics1.CumulativeReturn, metrics1.AnnualizedReturn, metrics1.SharpeRatio, 100/metrics1.MaxDrawdown];
    data2 = [metrics2.CumulativeReturn, metrics2.AnnualizedReturn, metrics2.SharpeRatio, 100/metrics2.MaxDrawdown];
    
    data = [data1; data2];
    
    data(isinf(data) | isnan(data)) = 0;
    
    if any(std(data) == 0)
       z_data = data;
    else
       z_data = zscore(data);
    end
    
    hf = figure('Name', '策略性能雷达图', 'Position', [200, 200, 700, 600]);
    lim_min = min(z_data(:)) - 0.1;
    lim_max = max(z_data(:)) + 0.1;
    if lim_min >= lim_max, lim_max = lim_min + 1; end
    create_spider(z_data', labels, [lim_min, lim_max], 5);
    
    title('策略性能雷达图 (Z-Score 标准化)');
    legend({metrics1.Name, metrics2.Name}, 'Location', 'southoutside');
    
    saveas(hf, fullfile(save_path, 'P2_Radar_Plot.emf'));
    fprintf('图2: Performance Attribution 雷达图 -> 已保存\n');
    
    function create_spider(P, fnames, lims, Ncirc)
        theta = linspace(0, 2*pi, size(P,1)+1)';
        P = [P; P(1,:)];
        polarplot(theta, ones(size(theta))*lims(2), 'k:'); hold on;
        for i = 1:Ncirc
            polarplot(theta, ones(size(theta))*(lims(1)+(lims(2)-lims(1))/Ncirc*i),'k:');
        end
        h = polarplot(theta, max(lims(1),P));
        rlim(lims);
        thetaticks(rad2deg(theta(1:end-1)));
        thetaticklabels(fnames);
    end
end

function plot_cumulative_alpha(T_Perf1, T_Perf2, save_path)
    [C, ia, ib] = intersect(T_Perf1.Time, T_Perf2.Time);
    if isempty(C), fprintf('图3: 策略Alpha累积图 -> 跳过 (无重叠日期)\n'); return; end
    
    T1_aligned = T_Perf1(ia, :);
    T2_aligned = T_Perf2(ib, :);
    
    alpha_ret = T1_aligned.R - T2_aligned.R;
    cum_alpha = cumsum(alpha_ret);
    
    hf = figure('Name', '累计Alpha对比图', 'Position', [300, 300, 1000, 500]);
    plot(C, cum_alpha, 'LineWidth', 2, 'DisplayName', '累计Alpha');
    hold on;
    
    slope = diff(cum_alpha);
    if length(slope) > 1
        [~, inflection_idx] = max(abs(diff(slope)));
        inflection_time = C(inflection_idx+1);
        inflection_val = cum_alpha(inflection_idx+1);
        scatter(inflection_time, inflection_val, 100, 'r', 'p', 'filled', 'DisplayName', 'Alpha拐点');
        text(inflection_time, inflection_val, ['  ' datestr(inflection_time, 'yyyy-mm-dd')], 'VerticalAlignment','bottom');
    end
    
    title(sprintf('累计Alpha (%s vs %s)', T_Perf1.Properties.Description, T_Perf2.Properties.Description));
    xlabel('时间');
    ylabel('累计超额收益 (Alpha)');
    legend('show', 'Location', 'northwest');
    grid on;
    datetick('x', 'yyyy-mm', 'keeplimits');
    
    saveas(hf, fullfile(save_path, 'P3_Cumulative_Alpha.emf'));
    fprintf('图3: 策略Alpha累积图 -> 已保存\n');
end

function plot_trade_density(TradeLog, TAsset, save_path)
    if isempty(TradeLog)
        strat_name = TradeLog.Properties.Description;
        fprintf('图4: 交易行为密度图 -> 跳过 (因为策略“%s”没有交易记录)\n', strat_name);
        return;
    end
    strat_name = TradeLog.Strategy{1};
    
    TAsset.DailyRange = TAsset.High - TAsset.Low;
    T_Vola = groupsummary(TAsset, 'Time', 'mean', 'DailyRange');
    T_Vola = renamevars(T_Vola, 'mean_DailyRange', 'Volatility');
    
    trade_dates = TradeLog.EntryDate;
    [is_member, loc] = ismember(trade_dates, T_Vola.Time);
    
    trade_vola = T_Vola.Volatility(loc(is_member));
    trade_dates_valid = trade_dates(is_member);

    hf = figure('Name', ['交易密度图: ' strat_name], 'Position', [400, 400, 1000, 600]);
    
    if isempty(trade_dates_valid)
        title(['交易密度图: ' strat_name ' (无有效交易可绘制)']);
        saveas(hf, fullfile(save_path, ['P4_Trade_Density_' strat_name '.emf']));
        fprintf('图4: 交易行为密度图 (%s) -> 已保存 (无有效交易)\n', strat_name);
        close(hf);
        return;
    end
    
    time_edges = linspace(min(T_Vola.Time), max(T_Vola.Time), 50);
    vola_min = min(T_Vola.Volatility);
    vola_max = max(T_Vola.Volatility);
    if vola_min == vola_max, vola_max = vola_min + 1; end
    vola_edges = linspace(vola_min, vola_max, 20);
    
    counts = histcounts2(datenum(trade_dates_valid), trade_vola, datenum(time_edges), vola_edges);
    
    imagesc(time_edges, vola_edges, counts');
    
    colorbar;
    ax = gca;
    ax.YDir = 'normal';
    
    title(['交易时机与频率图: ' strat_name]);
    xlabel('时间');
    ylabel('每日价格波动率');
    
    saveas(hf, fullfile(save_path, ['P4_Trade_Density_' strat_name '.emf']));
    fprintf('图4: 交易行为密度图 (%s) -> 已保存\n', strat_name);
end

function plot_rolling_efficiency(TradeLog1, TradeLog2, save_path, window_size)
    if isempty(TradeLog1) || isempty(TradeLog2)
        fprintf('图5: 滚动交易效率图 -> 跳过 (一个或多个策略没有交易记录)\n');
        return;
    end
    
    [eff1, time1, ~] = calc_rolling_eff(TradeLog1, window_size);
    [eff2, time2, ~] = calc_rolling_eff(TradeLog2, window_size);

    hf = figure('Name', '滚动交易效率图', 'Position', [500, 500, 1000, 600]);
    hold on;
    
    h_legend = [];

    if ~isempty(time1)
        h1 = plot(time1, eff1, '-b', 'LineWidth', 2, 'DisplayName', TradeLog1.Strategy{1});
        h_legend = [h_legend, h1];
    end
    if ~isempty(time2)
        h2 = plot(time2, eff2, '-g', 'LineWidth', 2, 'DisplayName', TradeLog2.Strategy{1});
        h_legend = [h_legend, h2];
    end
    
    if ~isempty(time1) && ~isempty(time2)
        [C, ia, ib] = intersect(time1, time2);
        if ~isempty(C)
            fill([C; flipud(C)], [eff2(ib); flipud(eff1(ia))], 'r', ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '效率差值');
        end
    end
    
    if ~isempty(h_legend), legend(h_legend, 'Location', 'northwest'); end
    
    hold off;
    
    title(['滚动交易效率 (窗口 = ' num2str(window_size) '笔交易)']);
    xlabel('时间');
    ylabel('滚动平均每笔交易盈亏');
    grid on;
    datetick('x', 'yyyy-mm', 'keeplimits');
    
    saveas(hf, fullfile(save_path, 'P5_Rolling_Efficiency.emf'));
    fprintf('图5: 滚动交易效率图 -> 已保存\n');

    function [eff, time, eff_std] = calc_rolling_eff(TLog, win)
        if height(TLog) < 2
            eff = []; time = []; eff_std = [];
            return;
        end
        TLog = sortrows(TLog, 'ExitDate');
        pnl = TLog.PnL_per_Unit;
        if length(pnl) < win, win = length(pnl); end
        eff = movmean(pnl, win, 'omitnan');
        eff_std = movstd(pnl, win, 'omitnan');
        time = TLog.ExitDate;
        
        valid_idx = ~isnan(eff);
        eff = eff(valid_idx); eff_std = eff_std(valid_idx); time = time(valid_idx);
    end
end