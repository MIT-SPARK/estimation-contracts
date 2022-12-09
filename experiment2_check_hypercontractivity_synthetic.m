%% Test Certifiable Hypercontractivity in robust rotation search (Wahba Problem)
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
Nset = [5 10 20 40 60 80 100 200]; % number of measurements % [5 10 20 40 60 80 100 200]
outrate = 0.5; % outlier rate (irrelevant in the computation of the hypercontractivity
Ct_pow_t_set = [1 2 3 4 5 6];
nrTests = 100; % number of runs for each configuration

%% allocate matrix to store results
results = zeros(length(Nset),length(Ct_pow_t_set)); 
runtime = zeros(length(Nset),length(Ct_pow_t_set)); 

%% generate random problems and test certifiable hyper-contractivity
for ind_N = 1:length(Nset)
    N = Nset(ind_N); % number of measurements
    for ind_Ct_set = 1:length(Ct_pow_t_set) 
        Ct_pow_t_input = Ct_pow_t_set(ind_Ct_set); % Ct for hypercontractivity
        nrHyperContractive = 0; % count number of certifiable hypercontractive instances
        totalTime = 0; % to measure time required by hyper contractivity check
        for test = 1:nrTests
            fprintf('=== test %d of %d, with N=%d, and outlier rate = %g ===\n',test,nrTests,N,outrate)
            [problem, R_gt] = createWahbaProblemData(N,outrate);
            [isHyperContractive,Q,Z,timeTest] = testHypercontractivity(problem.At,Ct_pow_t_input);
            totalTime = totalTime + timeTest;
            if isHyperContractive
                nrHyperContractive = nrHyperContractive + 1;
                fprintf('=== %d certifiable hypercontractive out of %d tests ===\n',nrHyperContractive,test)
            end
        end
        results(ind_N,ind_Ct_set) = nrHyperContractive;
        runtime(ind_N,ind_Ct_set) = totalTime/nrTests; % we record average time
    end
end

%save log_results_experiment2

%% plot/display results
figure
h = heatmap(Ct_pow_t_set,Nset,results);
% h.Title = 'Hypercontractive instances';
h.XLabel = '$C(t)^t$';
h.YLabel = 'Nr. measurements ($n$)';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
h.FontSize = 18;

%% timing figure
figure
h = heatmap(Ct_pow_t_set,Nset,runtime);
%h.Title = 'Average runtime';
h.XLabel = '$C(t)^t$';
h.YLabel = 'Nr. measurements ($n$)';
h.CellLabelFormat = '%.2f';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
h.FontSize = 18;








