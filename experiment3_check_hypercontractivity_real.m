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

Ct_pow_t_set = [1 2 3 4 5 6 7 8 9 10 15 20];

isHyperContractive_results = zeros(length(dataFileNames),length(Ct_pow_t_set));
timeTest_results = zeros(length(dataFileNames),length(Ct_pow_t_set));
eta_results = zeros(length(dataFileNames),1);

for fileInd=length(dataFileNames)
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
    plotImageStitchingResults
    
    %% store tightness
    eta_results(fileInd) = eta;
    
    for j = 1:length(Ct_pow_t_set)
        %% test hypercontractivity
        Ct_pow_t_input = Ct_pow_t_set(j); 
        [isHyperContractive,Q,Z,timeTest] = testHypercontractivity(problem.At,Ct_pow_t_input);
        isHyperContractive_results(fileInd,j) = isHyperContractive;
        timeTest_results(fileInd,j) = timeTest;
    end
    isHyperContractive_results
    timeTest_results
    eta_results
end

% save log_results_experiment3

%% print results for latex table
percentageHyperTest = 100*sum(isHyperContractive_results,1)/length(dataFileNames);
for i=1:length(percentageHyperTest)
    fprintf(' & $%.1f$ ', percentageHyperTest(i))
end
fprintf(' \n')





