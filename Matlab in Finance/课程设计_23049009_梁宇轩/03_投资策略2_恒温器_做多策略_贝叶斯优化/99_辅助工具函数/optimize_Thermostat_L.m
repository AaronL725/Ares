% AIM: 使用贝叶斯优化寻找恒温器策略的最优参数

%% P0, 清理与声明
clc; clear; close all;
global dataMap assetList params_fixed; % 定义全局变量

%% P1, 路径与数据准备

% 获取此脚本所在的文件夹
scriptFolder = fileparts(which(mfilename));
% 从脚本文件夹向上一级，找到项目根目录
projectRoot = fileparts(scriptFolder); 
% 将项目根目录及其所有子文件夹都添加到Matlab搜索路径
addpath(genpath(projectRoot));
cd(projectRoot); % 将当前工作目录切换到项目根目录
% ----------------------------------------------------

% 加载数据到全局变量
dataPath = fullfile(projectRoot, '01_数据及其函数');
dataMap = Data_Import_Thermostat(dataPath);
assetList = keys(dataMap);

% 过滤数据
startDate = datetime('2018-01-01');
endDate = datetime('2025-01-01');
for i = 1:length(assetList)
    assetTT = dataMap(assetList{i});
    dataMap(assetList{i}) = assetTT(timerange(startDate, endDate, 'closed'), :);
end
params_fixed.Lots = 1; % 定义固定的参数

%% P2, 定义优化变量与范围
optimVars = [
    optimizableVariable('swingTrendSwitch', [10, 50], 'Type', 'integer')
    optimizableVariable('swingPrcnt1', [0.2, 1.0], 'Type', 'real')
    optimizableVariable('swingPrcnt2', [0.2, 1.5], 'Type', 'real')
    optimizableVariable('atrLength', [5, 30], 'Type', 'integer')
    optimizableVariable('bollingerLengths', [20, 100], 'Type', 'integer')
    optimizableVariable('numStdDevs', [1.5, 3.5], 'Type', 'real')
    optimizableVariable('trendLiqLength', [20, 100], 'Type', 'integer')
];

%% P3, 运行贝叶斯优化
% 将目标函数句柄传递给 bayesopt
% MaxObjectiveEvaluations 设置了总的迭代次数，可以根据需要调整
results = bayesopt(@objective_Thermostat, optimVars, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'MaxObjectiveEvaluations', 50);

%% P4, 展示优化结果
bestParams = results.XAtMinObjective;
minNegativeSharpe = results.MinObjective;

fprintf('\n--- 贝叶斯优化完成 ---\n');
fprintf('找到的最佳夏普比率: %.4f\n', -minNegativeSharpe);
disp('对应的最佳参数组合为:');
disp(bestParams);