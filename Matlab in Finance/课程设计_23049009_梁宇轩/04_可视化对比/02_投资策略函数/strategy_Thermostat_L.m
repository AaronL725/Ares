% ===================== 1. 函数签名修改 =====================
% 增加了 params 输入参数，以接收外部传入的策略参数。
function [TPortfolio] = strategy_Thermostat_L(THistory, TCurrent, params)
% =======================================================

%% 1. 初始化投资组合
% ===================== 2. 输出表结构修改 =====================
% 在输出的 TPortfolio 表中增加了 'Time' 列，以符合框架要求。
TPortfolio = TCurrent(:, {'Time', 'AssetID'});
% =========================================================
TPortfolio.Weight = zeros(height(TCurrent), 1);
assetIDs = TCurrent.AssetID;

%% 2. 循环计算每个期货品种的交易信号
for i = 1:length(assetIDs)
    currentAssetID = assetIDs{i}; 
    TAssetHistory = THistory(strcmp(THistory.AssetID, currentAssetID), :);

    % ===================== 3. 参数获取方式修改 =====================
    % 从硬编码的固定值改为从 params 结构体中动态获取。
    swingTrendSwitch = params.swingTrendSwitch;
    swingPrcnt1 = params.swingPrcnt1;
    bollingerLengths = params.bollingerLengths;
    numStdDevs = params.numStdDevs;
    cmi_period = params.cmi_period;
    atrLength = params.atrLength;
    % ===========================================================
    
    if height(TAssetHistory) < max([cmi_period, bollingerLengths])
        continue;
    end
    
    %% 4. 计算指标 (逻辑保持不变)
    % (A) 潮汐指数 (CMI)
    cmi_close_diff = abs(TAssetHistory.Close(end) - TAssetHistory.Close(end - cmi_period + 1));
    cmi_price_range = max(TAssetHistory.High(end - cmi_period + 1 : end)) - min(TAssetHistory.Low(end - cmi_period + 1 : end));

    if cmi_price_range == 0
        cmiVal = 0;
    else
        cmiVal = cmi_close_diff / cmi_price_range * 100;
    end
    
    % (B) ATR for Swing Market
    tr = max(TAssetHistory.High - TAssetHistory.Low, max(abs(TAssetHistory.High - TAssetHistory.PrevClose), abs(TAssetHistory.Low - TAssetHistory.PrevClose)));
    myATR = mean(tr(end - atrLength + 1 : end));

    % (C) Bollinger Bands for Trend Market
    midLine = mean(TAssetHistory.Close(end - bollingerLengths + 1 : end));
    band = std(TAssetHistory.Close(end - bollingerLengths + 1 : end));
    dnBand = midLine - numStdDevs * band;

    %% 5. 生成交易决策 (逻辑保持不变)
    Position = 0; % 默认不持仓

    if cmiVal < swingTrendSwitch  % 震荡市
        % 当原始ATR突破信号触发时
        swingSellPt = TAssetHistory.Open(end) - swingPrcnt1 * myATR;
        if TAssetHistory.Close(end) < swingSellPt
            Position = 1;
        end
        
    else % 趋势市
        % 当布林通道突破信号触发时
        if TAssetHistory.Close(end) < dnBand
            Position = 1;
        end
    end
    
    %% 6. 分配权重
    TPortfolio.Weight(i) = Position;
end

end