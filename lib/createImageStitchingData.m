function [problem, R_gt, theta_gt] = createImageStitchingData(dataFileName)

load(dataFileName)

N = 70; % we max the points out at 70, as in the QUASAR paper
problem.N = N; % total number of measurements

%% create problem data
for i=1:N
    ai = Kinv * [double(matchedPtsImg1.Location(i,:)) 1]';
    a(:,i) = ai / norm(ai);
    bi = Kinv * [double(matchedPtsImg2.Location(i,:)) 1]';
    b(:,i) = bi / norm(bi);
end

problem.a = a;
problem.b = b;

%% from QUASAR paper
R_gt = R_quasar;
theta_gt = thetas;

%% set inlier bound
problem.betasq = 1e-3;
fprintf('betasq = %g\n',problem.betasq)

%% rearrange in matrix form, as in eq (4) of the paper 
problem.At = zeros(3,9,N);
for i=1:N
    problem.At(:,:,i) = kron( problem.a(:,i)' , eye(3) );
end