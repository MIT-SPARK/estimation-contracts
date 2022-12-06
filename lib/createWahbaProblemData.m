function [problem, R_gt, R_gt_list] = createWahbaProblemData(N,outrate,isAdversarial)
% isAdversarial = 1: creates 2 objects and adds similar noise on both
% isAdversarial = 2: creates 2 objects and adds adversarial noise on the first one
problem.N               = N; % total number of measurements
problem.Covariance      = 1e-4*eye(3); % noise covariance for inliers
problem.v1_distribution = 'uniform'; % sample uniformly from unit sphere
problem.nrOutliers      = round(N*outrate); % nr of outliers
problem.boundNoise      = true; % bound noise on inliers 
problem.normalize       = false; % don't normalize measurements to have unit norm (follows eq. (4) in the paper)

if nargin < 3
   isAdversarial = 0; 
end

if isAdversarial == 1
    problem.nrObjects = 2;
    problem.objectSize = [1-outrate,outrate];
    warning('Adversarial examples with consistent outliers')
end

if isAdversarial == 2
    problem.nrObjects = 2;
    problem.objectSize = [1-outrate,outrate];
    problem.Covariance = 1e-8*eye(3); % we make noise tiny, but add it back adversarially later on
    warning('Adversarial examples with consistent outliers and worst case noise')
end

if isAdversarial == 3
    problem.nrObjects = 5;
    problem.objectSize = 0.2 * ones(1,problem.nrObjects);
    warning('Adversarial examples with consistent outliers corresponding to 5 objects')
end

if isAdversarial > 3
   error('invalid choice of isAdversarial parameter') 
end

%% create problem data
[a, b, R_gt, problem]   = createWahbaProblem(problem);

%% potentially add back adversarial inlier noise
if isAdversarial == 2
    problem.Covariance = 1e-4*eye(3); % noise covariance for inliers
    problem.noiseBound = sqrt(problem.Covariance(2,2)*chi2inv(0.9999,3));
    nrInliers = N - problem.nrOutliers;    
    for k=1:nrInliers
        % add noise to inliers having worst case norm noiseBound (0.99 of max norm)
        noiseInliers = mvnrnd(zeros(1,3), problem.Covariance)'; 
        a(:,k) = a(:,k) + (0.99 * problem.noiseBound) * noiseInliers / norm(noiseInliers);
    end
end

%% set inlier bound
problem.betasq = problem.noiseBound^2
problem.a = a;
problem.b = b;

%% rearrange in matrix form, as in eq (4) of the paper 
problem.At = zeros(3,9,N);
for i=1:N
    problem.At(:,:,i) = kron( problem.a(:,i)' , eye(3) );
end

%% in the case of multiple objects, we only store 1
R_gt_list = R_gt;
if size(R_gt,3) > 1
    R_gt = R_gt(:,:,1); 
end 