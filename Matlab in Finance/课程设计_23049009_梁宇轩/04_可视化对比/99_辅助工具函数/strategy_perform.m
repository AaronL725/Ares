function T_Perform = strategy_perform(TPortfolio, TCurrent)

%% 1. 合并投资权重和行情数据
% 使用复合键 {'Time', 'AssetID'} 进行精确合并，避免内存问题。
% 这样可以确保某一天某一资产的权重，只匹配到同一天同一资产的行情。
T_Joined = outerjoin(TCurrent, TPortfolio, 'Keys', {'Time', 'AssetID'}, 'MergeKeys', true, 'Type', 'left');
T_Joined.Weight(isnan(T_Joined.Weight)) = 0; % 将未持仓的权重设为0

%% 2. 计算投资组合的日度收益率
% 这里使用 T_Joined 替代了原来的 TCurrent
dates = unique(T_Joined.Time);
daily_returns = zeros(length(dates), 1);

for i = 1:length(dates)
    TDay = T_Joined(T_Joined.Time == dates(i), :);
    
    % 获取当天持仓的资产
    TDay_Holding = TDay(TDay.Weight ~= 0, :);
    
    if isempty(TDay_Holding)
        daily_returns(i) = 0;
        continue;
    end

    % (A) 计算每个持仓合约的【当日盈亏金额 (P&L)】
    pnl_per_contract = (TDay_Holding.Close - TDay_Holding.PrevClose) .* TDay_Holding.Multi .* TDay_Holding.Weight;
    total_pnl = sum(pnl_per_contract, 'omitnan');

    % (B) 计算当日总保证金 (简化模型，使用前一日总名义价值作为分母)
    % 注：这是一个简化模型。真实保证金计算更复杂。
    total_nominal_value = sum(abs(TDay_Holding.PrevClose .* TDay_Holding.Multi), 'omitnan');
    
    if total_nominal_value == 0
        daily_returns(i) = 0;
    else
        daily_returns(i) = total_pnl / total_nominal_value;
    end
end

%% 3. 整理输出表格
T_Perform = table(dates, daily_returns, 'VariableNames', {'Time', 'R'});
T_Perform.RCum = cumprod(1 + T_Perform.R) - 1;

% 从 TPortfolio 的 UserData 中继承策略名称
if isprop(TPortfolio, 'Properties') && isfield(TPortfolio.Properties, 'Description') && ~isempty(TPortfolio.Properties.Description)
    T_Perform.Properties.Description = TPortfolio.Properties.Description;
end

end