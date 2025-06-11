function negativeSharpe = objective_Thermostat(paramsToOptimize)

    % 从全局变量获取数据，避免重复读取
    global TAsset_global;

    % 组合固定参数与待优化参数
    params = paramsToOptimize;
    % 在您的原始代码中，有一个名为 cmi_period 的参数被策略使用，但在优化变量中并未定义。
    % 这里我们暂时将其作为一个固定参数。如果它也需要优化，应将其添加到 optimizableVariable 列表中。
    params.cmi_period = 30; % 假设 cmi_period 是一个固定参数

    % --- 内部封装了完整的日度回测逻辑 ---
    all_trading_days = unique(TAsset_global.Time);
    nPeriod = length(all_trading_days);
    PeriodBegins = all_trading_days;
    history_window_size = 60;
    
    cellT_Perform = cell(nPeriod, 1);
    
    for iP = (history_window_size + 1):nPeriod
        current_day = PeriodBegins(iP);
        history_start_day = PeriodBegins(iP - history_window_size);
        history_end_day   = PeriodBegins(iP - 1);
        idx_History = (TAsset_global.Time >= history_start_day) & (TAsset_global.Time <= history_end_day);
        T_History = TAsset_global(idx_History, :);
        
        idx_Current = TAsset_global.Time == current_day;
        T_Current = TAsset_global(idx_Current, :);
        
        if isempty(T_Current) || isempty(T_History), continue; end
        
        % 调用【参数化】的策略函数
        TPortfolio = strategy_Thermostat_L(T_History, T_Current, params);
        
        % 度量策略表现 - 调用我们改进后的 strategy_perform.m 文件
        T_Perform = strategy_perform(TPortfolio, T_Current);
        
        cellT_Perform{iP} = T_Perform;
    end
    
    T_PerformAll = vertcat(cellT_Perform{:});
    
    % --- 计算最终评价指标 ---
    if isempty(T_PerformAll) || height(T_PerformAll) < 20 % 如果交易日太少，则认为无效
        negativeSharpe = 100; % 返回一个很大的“惩罚”值
        return;
    end
    
    T_PerformAll.RCum = cumprod(1 + T_PerformAll.R);
    
    annualReturn = T_PerformAll.RCum(end)^(252/height(T_PerformAll)) - 1;
    annualVolatility = std(T_PerformAll.R) * sqrt(252);
    
    if annualVolatility < 1e-6 % 如果波动为0，避免除零错误
        sharpeRatio = 0;
    else
        sharpeRatio = annualReturn / annualVolatility;
    end
    
    % 贝叶斯优化默认寻找最小值，因此返回负的夏普比率
    negativeSharpe = -sharpeRatio;

end