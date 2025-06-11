function [T_Asset, T_SH_Index] = Data_Import_Thermostat(dataPath)

%% 1. 读取并处理期货数据
futureFile = fullfile(dataPath, 'FutureData_23049009_梁宇轩_基于市场状态自适应的恒温器期货量化交易策略与贝叶斯参数优化研究.csv');
opts = detectImportOptions(futureFile);
opts.VariableTypes{strcmp(opts.VariableNames, 'Asset_ID')} = 'char';
opts.VariableTypes{strcmp(opts.VariableNames, 'Datetime')} = 'datetime';
opts = setvaropts(opts, 'Datetime', 'InputFormat', 'yyyy-MM-dd');
T_Asset = readtable(futureFile, opts);
T_Asset = renamevars(T_Asset, ["Asset_ID", "Datetime"], ["AssetID", "Time"]);

%% 2. 计算附加列
T_Asset = sortrows(T_Asset, {'AssetID', 'Time'});
assetGroups = findgroups(T_Asset.AssetID);

% 计算前一日收盘价 PrevClose
prevCloses = splitapply(@(x) {[NaN; x(1:end-1)]}, T_Asset.Close, assetGroups);
T_Asset.PrevClose = vertcat(prevCloses{:});

% 计算简单收益率 R = (Close / PrevClose) - 1
T_Asset.R = (T_Asset.Close ./ T_Asset.PrevClose) - 1;

% 计算平均价格 P
T_Asset.P = (T_Asset.Open + T_Asset.High + T_Asset.Low + T_Asset.Close) / 4;

%% 3. 清理数据
T_Asset(isnan(T_Asset.PrevClose) | isinf(T_Asset.R) | isnan(T_Asset.R), :) = [];

%% 4. 读取上证指数数据
sh_index_file = fullfile(dataPath, 'SH_Index_23049009_梁宇轩_基于市场状态自适应的恒温器期货量化交易策略与贝叶斯参数优化研究.csv');
opts_sh = detectImportOptions(sh_index_file);
opts_sh = setvartype(opts_sh, 'Datetime', 'datetime');
opts_sh = setvaropts(opts_sh, 'Datetime', 'InputFormat', 'yyyy-MM-dd');
T_SH_Index = readtable(sh_index_file, opts_sh);
T_SH_Index = renamevars(T_SH_Index, "Datetime", "Time");

end