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
sigAB = sigA*sigB*rhoAB;

wA = 0:0.1:1;
wB = 1-wA; 

% calculation
muP = wA*muA + wB*muB;
sigP = sqrt(wA.^2*(sigA^2) + wB.^2*(sigB^2) + 2.*wA.*wB*sigAB);

wAmin = (sigB^2 - sigAB)/(sigA^2 + sigB^2 - 2*sigAB);
wBmin = 1 - wAmin;
muMin = wAmin*muA + wBmin*muB;
sigMin = sqrt(wAmin.^2*(sigA^2) + wBmin.^2*(sigB^2) + 2.*wAmin.*wBmin*sigAB);

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

title('Minimal Risk Portfolio: $\sigma_{min}$', 'Interpreter', 'latex');
set(0, 'DefaultAxesFontName', 'Times New Roman');
xlabel('$\sigma_P$', 'Interpreter', 'latex');
ylabel('Mean $\mu_P$', 'Interpreter', 'latex');
