function [TPortfolio] = strategy_Thermostat_L(THistory, TCurrent, params)

%% 1. 初始化投资组合
TPortfolio = TCurrent(:, {'AssetID'});
TPortfolio.Weight = zeros(height(TCurrent), 1);
assetIDs = TCurrent.AssetID;

%% 2. 循环计算每个期货品种的交易信号
for i = 1:length(assetIDs)
    currentAssetID = assetIDs{i}; 
    TAssetHistory = THistory(strcmp(THistory.AssetID, currentAssetID), :);

    %%for i = minLength+1:n 3. 从 params 结构体中解构策略参数
    swingTrendSwitch = params.swingTrendSwitch;
    swingPrcnt1 = params.swingPrcnt1;
    bollingerLengths = params.bollingerLengths;
    numStdDevs = params.numStdDevs;
    cmi_period = params.cmi_period;
    atrLength = params.atrLength;
    
    if height(TAssetHistory) < max([cmi_period, bollingerLengths, atrLength])
        continue;
    end
    
    %% 4. 计算指标
    cmi_close_diff = abs(TAssetHistory.Close(end) - TAssetHistory.Close(end - cmi_period + 1));
    cmi_price_range = max(TAssetHistory.High(end - cmi_period + 1 : end)) - min(TAssetHistory.Low(end - cmi_period + 1 : end));
    if cmi_price_range == 0, cmiVal = 0; else, cmiVal = cmi_close_diff / cmi_price_range * 100; end
    
    tr = max(TAssetHistory.High - TAssetHistory.Low, max(abs(TAssetHistory.High - TAssetHistory.PrevClose), abs(TAssetHistory.Low - TAssetHistory.PrevClose)));
    myATR = mean(tr(end - atrLength + 1 : end));

    midLine = mean(TAssetHistory.Close(end - bollingerLengths + 1 : end));
    band = std(TAssetHistory.Close(end - bollingerLengths + 1 : end));
    dnBand = midLine - numStdDevs * band;

    %% 5. 生成交易决策
    Position = 0;
    if cmiVal < swingTrendSwitch
        swingSellPt = TAssetHistory.Open(end) - swingPrcnt1 * myATR;
        if TAssetHistory.Close(end) < swingSellPt
            Position = 1;
        end
    else
        if TAssetHistory.Close(end) < dnBand
            Position = 1;
        end
    end
    
    %% 6. 更新投资组合权重
    idx = strcmp(TPortfolio.AssetID, currentAssetID);
    TPortfolio.Weight(idx) = Position;
end

end