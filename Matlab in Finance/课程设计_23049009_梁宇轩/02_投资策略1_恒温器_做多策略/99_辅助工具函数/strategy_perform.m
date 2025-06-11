function T_Perform = strategy_perform(TPortfolio, TCurrent)

%% 1. 合并投资权重和行情数据
TCurrent = outerjoin(TCurrent, TPortfolio, 'Keys', 'AssetID', 'MergeKeys', true, 'Type', 'left');
TCurrent.Weight(isnan(TCurrent.Weight)) = 0;

%% 2. 计算投资组合的日度收益率
dates = unique(TCurrent.Time);
daily_returns = zeros(length(dates), 1);

for i = 1:length(dates)
    TDay = TCurrent(TCurrent.Time == dates(i), :);
    
    % 获取当天持仓的资产
    TDay_Holding = TDay(TDay.Weight ~= 0, :);
    
    if isempty(TDay_Holding)
        daily_returns(i) = 0;
        continue;
    end

    % (A) 计算每个持仓合约的【当日盈亏金额 (P&L)】
    % 通过乘以Weight，自动处理多空方向：
    % Weight=1 (做多): (Close-PrevClose)*Multi
    % Weight=-1(做空): (Close-PrevClose)*Multi*(-1) = (PrevClose-Close)*Multi
    pnl_per_contract = (TDay_Holding.Close - TDay_Holding.PrevClose) .* TDay_Holding.Multi .* TDay_Holding.Weight;
    total_pnl = sum(pnl_per_contract, 'omitnan');

    % (B) 计算投资组合在【期初】的【总名义价值】
    total_notional_value = sum(TDay_Holding.PrevClose .* TDay_Holding.Multi .* abs(TDay_Holding.Weight), 'omitnan');
    
    % (C) 计算组合的【日收益率】
    if total_notional_value == 0
        daily_returns(i) = 0;
    else
        daily_returns(i) = total_pnl / total_notional_value;
    end
end

%% 3. 格式化输出
T_Perform = table;
T_Perform.Time = TCurrent.Time(1);
T_Perform.R = daily_returns(1);

end