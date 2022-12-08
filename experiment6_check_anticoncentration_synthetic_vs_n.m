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
Nset = [10 50 100 200]; % number of measurements
outrateSet = [0 0.1 0.3 0.5 0.7 0.9]; % outlier rate
nrTests = 10; % number of runs for each configuration
isModifiedWahba = 0;

alphaBar = 0.8;
eta = (2/alphaBar) * (2 - 2* (1-alphaBar)/(alphaBar)) % smallest eta such that bound is informative

%% matrix storing results
results = zeros(length(Nset),length(outrateSet), nrTests); 
runtime1 = zeros(length(Nset),length(outrateSet), nrTests); 
runtime2 = zeros(length(Nset),length(outrateSet), nrTests); 

%% generate random problems and test anti-concentration
for N_ind = 1:length(Nset)
    N = Nset(N_ind); % number of measurements
    for outrate_ind = 1:length(outrateSet)
        outrate = outrateSet(outrate_ind); % outlier rate
        nrAntiConcentrated = 0;
        for test = 1:nrTests
            fprintf('=== test %d of %d, with N=%d, and outlier rate = %g ===\n',test,nrTests,N,outrate)
            [problem, R_gt] = createWahbaProblemData(N,outrate);
            problem.outrate = outrate;
            [isAntiConcentrated,timeTest1,timeTest2] = testAntiConcentration(problem,eta,isModifiedWahba);
            
            if isAntiConcentrated
                nrAntiConcentrated = nrAntiConcentrated + 1;
            end
            
            results(N_ind,outrate_ind,test) = isAntiConcentrated;
            runtime1(N_ind,outrate_ind,test) = timeTest1;
            runtime2(N_ind,outrate_ind,test) = timeTest2;
            fprintf('================ %d sos out of %d ====================\n',nrAntiConcentrated,test)
        end
    end
end

% save log_results_experiment6_eta_at_080 % good results for isModified = 1

%% plot/display results
resultsPercentage = 100*sum(results,3)/nrTests;
aveRuntime1 = mean(runtime1,3);
aveRuntime2 = mean(runtime2,3);

aveAveRuntime1 = mean(runtime1(:)) % overall average
aveAveRuntime2 = mean(runtime2(:)) % overall average

figure
h = heatmap(outrateSet,Nset,resultsPercentage);
% h.Title = 'Hypercontractive instances';
h.XLabel = 'Outlier rate ($\beta$)';
h.YLabel = 'Nr. measurements ($n$)';
h.CellLabelFormat = '%.2f';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
h.FontSize = 18;

%% timing figure
figure
h = heatmap(outrateSet,Nset,aveRuntime1);
%h.Title = 'Average runtime';
h.XLabel = 'Outlier rate $\beta$';
h.YLabel = 'Nr. measurements ($n$)';
h.CellLabelFormat = '%.2f';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
h.FontSize = 18;

figure
h = heatmap(outrateSet,Nset,aveRuntime2);
%h.Title = 'Average runtime';
h.XLabel = 'Outlier rate $\beta$';
h.YLabel = 'Nr. measurements ($n$)';
h.CellLabelFormat = '%.2f';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
h.FontSize = 18;


