%% Test Certifiable Anti-Concentration in robust rotation search (Wahba Problem)
%% Synthetic data
%% Paper: "Estimation Contracts for Outlier-Robust Geometric Perception"
%% Luca Carlone, Nov 7, 2022
%% example data from: https://github.com/MIT-SPARK/CertifiablyRobustPerception/

clc
clear
close all
restoredefaultpath
addSpecificPaths

%% set parameters for Monte Carlo runs
N = 50;
alphaBarSet = [0.55 0.6 0.7 0.8]; % inlier rate for outrateSet = [0.2 0.3 0.4 0.45]; % outlier rate
for i=1:length(alphaBarSet)
    alphaBar = alphaBarSet(i);
    etaSet(i) = (2/alphaBar) * (2 - 2* (1-alphaBar)/(alphaBar)) % smallest eta such that bound is informative
end
outrateSet = [0 0.1 0.3 0.5 0.7 0.9]; % outlier rate
nrTests = 10; % number of runs for each configuration
isModifiedWahba = 1;

%% matrix storing results
results = zeros(length(etaSet),length(outrateSet), nrTests); 
runtime1 = zeros(length(etaSet),length(outrateSet), nrTests); 
runtime2 = zeros(length(etaSet),length(outrateSet), nrTests); 

%% generate random problems and test anti-concentration
for eta_ind = 1:length(etaSet)
    eta = etaSet(eta_ind); % number of measurements
    for outrate_ind = 1:length(outrateSet)
        outrate = outrateSet(outrate_ind); % outlier rate
        nrAntiConcentrated = 0;
        for test = 1:nrTests
            fprintf('=== test %d of %d, with eta=%g, and outlier rate = %g ===\n',test,nrTests,eta,outrate)
            [problem, R_gt] = createWahbaProblemData(N,outrate);
            problem.outrate = outrate;
            [isAntiConcentrated,timeTest1,timeTest2] = testAntiConcentration(problem,eta,isModifiedWahba);
                 
            % log results
            results(eta_ind,outrate_ind,test) = isAntiConcentrated;
            runtime1(eta_ind,outrate_ind,test) = timeTest1;
            runtime2(eta_ind,outrate_ind,test) = timeTest2;
            
            % display partial results
            if isAntiConcentrated
                nrAntiConcentrated = nrAntiConcentrated + 1;
            end
            fprintf('================ %d sos out of %d ====================\n',nrAntiConcentrated,test)
        end
    end
end

%% plot/display results
resultsPercentage = 100*sum(results,3)/nrTests;
aveRuntime1 = mean(runtime1,3);
aveRuntime2 = mean(runtime2,3);

aveAveRuntime1 = mean(runtime1(:)) % overall average
aveAveRuntime2 = mean(runtime2(:)) % overall average

save log_results_experiment7_isModified_1

figure
h = heatmap(outrateSet,etaSet,resultsPercentage);
% h.Title = 'Hypercontractive instances';
h.XLabel = 'Outlier rate ($\beta$)';
h.YLabel = '$\eta$';
h.YDisplayLabels = compose('%.2f',str2double(h.YDisplayLabels));
h.CellLabelFormat = '%.2f';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
h.FontSize = 18;

%% timing figure
figure
h = heatmap(outrateSet,etaSet,aveRuntime1);
%h.Title = 'Average runtime';
h.XLabel = 'Outlier rate $\beta$';
h.YLabel = '$\eta$';
h.YDisplayLabels = compose('%.2f',str2double(h.YDisplayLabels));
h.CellLabelFormat = '%.2f';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
h.FontSize = 18;

figure
h = heatmap(outrateSet,etaSet,aveRuntime2);
%h.Title = 'Average runtime';
h.XLabel = 'Outlier rate $\beta$';
h.YLabel = '$\eta$';
h.YDisplayLabels = compose('%.2f',str2double(h.YDisplayLabels));
h.CellLabelFormat = '%.2f';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
h.FontSize = 18;


