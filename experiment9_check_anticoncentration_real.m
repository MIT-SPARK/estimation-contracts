%% Test Certifiable Hypercontractivity in robust rotation search (Wahba Problem)
%% Real data from image stitching application
%% Paper: "Estimation Contracts for Outlier-Robust Geometric Perception"
%% Luca Carlone, Nov 7, 2022
%% example data provided by Heng Yang (same as the one used in the QUASAR paper: https://arxiv.org/pdf/1905.12536.pdf)

clc
clear
close all
restoredefaultpath
addSpecificPaths

%% name of mat file containing images and feature data
dataFileNames = {...
    'data/stitching_1_7.mat',...
    'data/stitching_7_13.mat',... %-x
    'data/stitching_13_19.mat',...
    'data/stitching_19_25.mat',...
    'data/stitching_25_31.mat',...
    'data/stitching_31_37.mat',...
    'data/stitching_37_43.mat',... %-x
    'data/stitching_43_49.mat',... %-x
    'data/stitching_49_55.mat',...
    'data/stitching_55_61.mat',...%-x
    'data/stitching_61_67.mat',...
    'data/stitching_67_1.mat'} %-xx

alphaBarSet = [0.55 0.6 0.7 0.8]; % inlier rate for outrateSet = [0.2 0.3 0.4 0.45]; % outlier rate
for i=1:length(alphaBarSet)
    alphaBar = alphaBarSet(i);
    etaSet(i) = (2/alphaBar) * (2 - 2* (1-alphaBar)/(alphaBar)) % smallest eta such that bound is informative
end

isAnticoncentrated_results = zeros(length(dataFileNames),length(etaSet));
timeTest1_results = zeros(length(dataFileNames),length(etaSet));
timeTest2_results = zeros(length(dataFileNames),length(etaSet));
eta_results = zeros(length(dataFileNames),1);

for fileInd=1:length(dataFileNames)
    dataFileName = dataFileNames{fileInd}
    
    %% parse images to get real data
    [problem, R_gt, theta_gt] = createImageStitchingData(dataFileName);
    
    %% contrarily to experiment1, where we know the ground truth inliers, here we estimate them
    [R_stride,theta_stride,eta] = solveWahba(problem.a,problem.b,problem.betasq,...
        R_gt,... % R_gt only used for performance evaluation
        sdpnalpath,manoptpath,stridepath); % paths to solvers
    
    %% build inliers/outliers (estimated from real data this time)
    problem.outliersIds = find(theta_stride < 0);

    %% check results against QUASAR results (marked as "gt")
    theta_diff = abs(theta_stride - theta_gt);
    diffOutliers = sum(theta_diff>0); % Sdp formulation is slightly different from QUASAR,..
    % ... but number of inliers should roughly match (+/-1)
    R_err = R_gt' * R_stride; 
    axang = rotm2axang(R_err);
    angleErrDeg = rad2deg(axang(4)); % small, as expected (seems rotm2axang approximates to zero)
    fprintf('diffOutliers = %d, rotError[deg] = %g\n',diffOutliers,angleErrDeg)
    
    %% visualize image stitching results
    % plotImageStitchingResults
    
    %% store tightness
    eta_results(fileInd) = eta;

    %% use outliers from stride
    problem.outliersIds = find(theta_stride < 0);
    problem.nrOutliers = length(problem.outliersIds);
    problem.outrate = problem.nrOutliers / problem.N
    
    for j = 1:length(etaSet)
        %% test hypercontractivity
        eta = etaSet(j); 
        [isAntiConcentrated,timeTest1,timeTest2] = testAntiConcentration(problem,eta);
        isAnticoncentrated_results(fileInd,j) = isAntiConcentrated;
        timeTest1_results(fileInd,j) = timeTest1;
        timeTest2_results(fileInd,j) = timeTest2;
    end
    isAnticoncentrated_results
    timeTest1_results
    timeTest2_results
    eta_results
end

% save log_results_experiment8

%% print results for latex table
percentageHyperTest = 100*sum(isAnticoncentrated_results,1)/length(dataFileNames);
for i=1:length(percentageHyperTest)
    fprintf(' & $%.1f$ ', percentageHyperTest(i))
end
fprintf(' \n')





