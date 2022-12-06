function [isAntiConcentrated,timeTest1,timeTest2] = testAntiConcentration(problem, eta_input,isModifiedWahba)

%% set (C, delta, M)-certifiable anti-concentration parameters according to Prop 12:
N = problem.N;
inliersIds = setdiff([1:N], problem.outliersIds); % anti-concentration only involves the inliers
nrInliers = N - problem.nrOutliers;
alpha = 1-problem.outrate; % inlier rate
barc = sqrt(problem.betasq); % maximum error for inliers

Mx = sqrt(3); % bound on ||x||
M = 2*Mx; % bound on || x - x_gt ||

if nargin < 2
    alphaBar = 0.7;
    eta = (2/alphaBar) * (2 - 2* (1-alphaBar)/(alphaBar)); % smallest eta such that bound is informative
else
    warning('using input eta')
    eta = eta_input
end
if nargin<3
    isModifiedWahba = 0;
end

delta = 2*barc;
C = alpha^2 * eta^2 * (1-2*barc)^2 / (32*barc); 
C_delta_Msq = C*delta*M^2

%% Matrices that have to satisfy property
At = problem.At;

if isModifiedWahba
    a = problem.a;
    At = zeros(9,9,N);
    for i=1:N
        a2 = (eye(3) - a(:,i) * a(:,i)') * rand(3,1);
        a2 = a2 / norm(a2);
        % check orthogonality and unit norm:
        if abs( norm(a2) - 1) > 1e-4 || abs( a(:,i)' * a2 ) > 1e-4
            i
            norm(a2)
            abs( a(:,i)' * a2 )
            error('error in a2 generation')
        end
        a3 = cross(a(:,i),a2);
        At(:,:,i) = [kron( a(:,i)' , eye(3) ) ;  kron( a2' , eye(3) )  ;  kron( a3' , eye(3) )];
        
        assert( norm(At(:,:,i) * At(:,:,i)' - eye(9)) <1e-6 ); % check that matrix is orthogonal
    end
end

%% coefficients of polynomial p (assumed degree 2): p(x) = 1 + c2*x^2
c2 = -0.1; 

%% prepare variables for sostools (to check sos proofs)
pvar v1 v2 v3 v4 v5 v6 v7 v8 v9
v = [v1;v2;v3;v4;v5;v6;v7;v8;v9];

%% check anti-concentration -- condition 1:
A = [];
degree = 4; % for this is typically enough to check deg 4
isAntiConcentrated = 1;
timeTest1 = 0;
nrRuns = 0;
for i=1:length(inliersIds)
    ind = inliersIds(i);
    normsq_Ai_v(i) = sum( (At(:,:,ind) * v ).^2 ); % ||Ai' v||^2 
    
    %% check anti-concentration -- condition 1:
    %p(i) = 1 + c2 * normsq_Ai_v(i) + c4 * normsq_Ai_v(i)^2; % c0 = 1 (by definition)
    p(i) = 1 + c2 * normsq_Ai_v(i);
    
    tstart = tic;
    [gam1,vars,opt] = findbound( p(i)^2 - (1-delta)^2, [delta^2-normsq_Ai_v(i)],[],degree);
    timeTest1 = timeTest1+toc(tstart);
    nrRuns = nrRuns + 1;
    if gam1 < 0
        isAntiConcentrated = 0;
        warning('found an Ai that is not anti-concentrated (condition 1')
    end
    gam1
    warning('remove i==1 if')
    
    A = [ A;At(:,:,ind) ];
end
timeTest1 = timeTest1/nrRuns;

if isAntiConcentrated == 1 % if we haven't found a case already making the condition false
    %% check anti-concentration -- condition 2:
    degree = 6; % here we need at least degree 6
    normsq_v = sum(v.^2);
    
    tstart = tic;
    [gam2,vars,opt] = findbound( C_delta_Msq - normsq_v * (1/nrInliers) * sum(p.^2),[M^2 - normsq_v],[],degree);
    timeTest2 = toc(tstart);
    gam2
    if gam2 < 0
        isAntiConcentrated = 0;
        warning('not anti-concentrated (condition 2')
    end
end

