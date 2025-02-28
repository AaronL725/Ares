% clear the workspace
clc;
close all;
clear; 

% input
muA = 0.1775;
muB = 0.055;
sigA = sqrt(0.067);
sigB = sqrt(0.013);
rhoAB = -0.164;
rf = 0.05;
sigAB = sigA*sigB*rhoAB;

wA = 0:0.1:1;
wB = 1-wA; 

% function
function [sigP, muP, wA, wB] = FrontierWeight(wA, muA, muB, sigA, sigB, sigAB)
    arguments (Input)
        wA (:,:) double
        muA (1,1) double
        muB (1,1) double
        sigA (1,1) double
        sigB (1,1) double
        sigAB (1,1) double
    end

    arguments (Output)
        sigP (:,1) double
        muP (:,1) double
        wA (:,1) double
        wB (:,1) double
    end

    % wA = 0:0.1:1;
    wB = 1 - wA;

    % calculation
    muP = wA*muA + wB*muB;
    sigP = sqrt(wA.^2*(sigA^2) + wB.^2*(sigB^2) + 2.*wA.*wB*sigAB);
end

% calculation
muP = wA*muA + wB*muB;
sigP = sqrt(wA.^2*(sigA^2) + wB.^2*(sigB^2) + 2.*wA.*wB*sigAB);

wAmin = (sigB^2 - sigAB)/(sigA^2 + sigB^2 - 2*sigAB);
wBmin = 1 - wAmin;
muMin = wAmin*muA + wBmin*muB;
sigMin = sqrt(wAmin.^2*(sigA^2) + wBmin.^2*(sigB^2) + 2.*wAmin.*wBmin*sigAB);

wAtanTop = (muA - rf)*sigB^2 - (muB - rf)*sigAB;
wAtanBot = (muB - rf)*sigA^2 + (muA - rf)*sigB^2 - (muA - rf + muB - rf)*sigAB;
wAtan = wAtanTop/wAtanBot;
[sigPtan, muPtan, wAtan, wBtan] = FrontierWeight(wAtan, muA, muB, sigA, sigB, sigAB);

% plot
figure(1)
hold on;

plot(sigP, muP, 'b.-');

plot(sigA, muA, 'rh', 'MarkerFaceColor', 'r');
text(sigA + 0.005, muA, '$A$', 'Interpreter', 'latex');
plot(sigB, muB, 'ro');
text(sigB + 0.005, muB, '$B$', 'Interpreter', 'latex');

plot(sigMin, muMin, 'ks', 'MarkerFaceColor', 'g');
text(sigMin + 0.005, muMin, '$P_{min}$', 'Interpreter', 'latex');

set(0, 'DefaultAxesFontName', 'Times New Roman');
xlabel('$\sigma_P$', 'Interpreter', 'latex');
ylabel('Mean $\mu_P$', 'Interpreter', 'latex');

plot(sigPtan, muPtan, 'ko', 'MarkerFaceColor', 'k');
plot(0, rf, 'ko', 'MarkerFaceColor', 'k');
plot([sigPtan, 0], [muPtan, rf], 'k-');

text(sigPtan + 0.005, muPtan, '$P_{tan}$', 'Interpreter', 'latex');
title('Tangence Portfolio', 'Interpreter', 'latex');
