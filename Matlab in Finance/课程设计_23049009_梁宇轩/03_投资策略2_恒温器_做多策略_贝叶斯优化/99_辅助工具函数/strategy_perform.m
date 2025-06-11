function T_Perform = strategy_perform(TPortfolio, TCurrent)
% STRATEGY_PERFORM 根据给定的投资组合和当日行情计算策略的单日表现。
%
% 该函数采用期货中常用的方法：计算当日总盈亏(P&L)，然后除以当日部署的总名义资本，
% 得到策略的日收益率。
%
% 输入:
%   - TPortfolio: 包含 {'AssetID', 'Weight'} 的投资组合表。Weight=1为做多，-1为做空。
%   - TCurrent:   包含单日市场行情的表，需含有 
%                 {'AssetID', 'Time', 'Close', 'PrevClose', 'Multi'} 列。
%
% 输出:
%   - T_Perform:  一个 1x2 的表，包含 {'Time', 'R'} 两列，R为策略当日收益率。

%% 1. 数据准备与输入验证
% 确保 TCurrent 表中有数据，否则返回零收益。
if height(TCurrent) < 1
    T_Perform = table(NaT, 0, 'VariableNames', {'Time', 'R'});
    return;
end
current_day = TCurrent.Time(1); % 获取当前日期

%% 2. 合并投资组合权重与市场行情数据
% 使用 outerjoin 将权重附加到行情数据上。
% 投资组合中未包含的资产，其权重为 NaN，我们将其置为 0。
TDay = outerjoin(TCurrent, TPortfolio, 'Keys', 'AssetID', 'MergeKeys', true, 'Type', 'left');
TDay.Weight(isnan(TDay.Weight)) = 0;

%% 3. 向量化计算策略表现 (无循环)

% (A) 计算每个合约的当日盈亏 (Profit & Loss)
% 公式: (当日收盘价 - 昨日收盘价) * 合约乘数 * 仓位权重
% TDay.Weight 自动处理了多空方向。
pnl_per_contract = (TDay.Close - TDay.PrevClose) .* TDay.Multi .* TDay.Weight;
total_pnl = sum(pnl_per_contract, 'omitnan');

% (B) 计算策略部署的总名义资本 (Total Capital Deployed)
% 公式: abs(仓位权重) * 当日收盘价 * 合约乘数
% 使用 abs() 是因为无论多空，占用的名义资本都是正值。
capital_per_contract = abs(TDay.Weight) .* TDay.Close .* TDay.Multi;
total_capital = sum(capital_per_contract, 'omitnan');

% (C) 计算当日组合收益率
if total_capital > 1e-9 % 使用一个极小值以避免浮点数除零错误
    daily_return = total_pnl / total_capital;
else
    % 如果没有部署资金 (即空仓), 当日收益为 0。
    daily_return = 0;
end

%% 4. 构建标准输出表
T_Perform = table(current_day, daily_return, 'VariableNames', {'Time', 'R'});

end