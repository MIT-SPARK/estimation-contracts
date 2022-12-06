%% solve SDP relaxation for List Decodable Regression using cvx (with w in {0;1})
% Algorithm 5 in the paper "Estimation Contracts for Outlier-Robust Geometric Perception"
% but using a sparse order 2 moment relaxation

v1=problem.a;
v2=problem.b;
betasq=problem.betasq;

N = size(v1,2); % nr vectors (measurements)
%% sparse monomial basis
% X = [1 w_1 ... w_N q\tran, q_1\tran, q_2\tran, ..., q_N\tran] has size n = (N+1)+4*(N+1)
n = (N+1)+4*(N+1);

%% solve SDP relaxation
% w_i || v1 - R*v1 ||^2 <= barc2 constraints
[A,b] = maxResidualConstraints(v1,v2,betasq);

% sum_i w_i == alpha * N
if recoverAllHypotheses == 0
    alpha = 1-outrate;
else
    alpha = min(outrate,1-outrate)
end
assert(abs(round(alpha*N) - alpha*N)<1e-7, 'alpha*N is not integer'); % check it's integer
[Aalpha,balpha] = fixedNumberInliers(N,alpha);

% check that constraints are correct using gt as a feasible solution
unitTest(R_gt,problem.theta_gt(:),N,alpha,v1,v2,betasq,A,b,Aalpha,balpha,isAdversarial);

offset = N+1; % shortcut to access part of the moment matrix
q_ind = [1:4]+offset;

tstart = tic;
    
%% solve relaxation
cvx_begin % sdp
cvx_solver mosek
variable X(n,n) symmetric
minimize( norm( X(2:N+1,1) ) ) % norm(w) % quad_form( X(2:N+1,1), eye(N) ) ) % norm(w)
subject to
X == semidefinite(n);
X(1,1) == 1;
trace(X(q_ind,q_ind)) == 1;
for k=1:N
    %% moment constraints
    % top left diagonal (w_i^2 = w_i)
    X(k+1,k+1) == X(1,k+1);
    % bottom right block diagonal (q_i q_i' = w_i^2 q q' = w_i q q' = q q_i)
    idx = blkIndices(k+1,4)+offset;
    X(idx,idx) == X(q_ind,idx);
    % intermediate entries (q_i' = w_i q')
    X(1,idx) == X(k+1,q_ind);
    % intermediate entries (w_i q_i' = w_i^2 q = w_i q = q_i)
    X(k+1,idx) == X(1,idx);
    
    %% max residual constraint
    trace( A(:,:,k) * X ) <= b(k); % max residual constraints 
end
%% structural constraints
for k1=1:N
    for k2=k1+1:N+1
        % symmetric off diagonal
        idx1 = blkIndices(k1,4)+offset;
        idx2 = blkIndices(k2,4)+offset;
        X(idx1,idx2) == X(idx1,idx2)';
        % tie theta_i * theta_j to trace of q_i * q_j' (including theta_i vs q * q_i
        X(k1,k2) == trace( X(idx1,idx2) );
    end
end
%% desired fraction of inliers
trace( Aalpha * X ) == balpha;
cvx_end

%% inspect optimal solution of SDP
f_sdp   = cvx_optval; % lower bound
e = eig(X);
disp('top 5 eigenvalues')
e(end-5:end)
if isnan(f_sdp)==1
    warning('cvx failed to solve correctly')
end

%% Compute relaxation gap
f_est   = sqrt(alpha*N); % optimal cost of original problem (with simple norm) is always norm of w and w has apha*N nonzero entries = 1 
subopt  = abs(f_est - f_sdp) / (1+abs(f_est)+abs(f_sdp));
eta_norm = subopt;

% optimal cost of original problem (with squared norm) is always alpha n 
eta  = abs(f_est^2 - f_sdp^2) / (1+abs(f_est^2)+abs(f_sdp^2));

%% Rounding (similar to Algorith 5 in the estimation contracts paper,
% but solution is also normalized here
PExp_w = X([2:N+1],1);
PExp_xw = zeros(4,N);
v = zeros(4,N); % vectors v as in Algorithm 5
p = zeros(N,1); % probabilities
for k=1:N
    idx = blkIndices(k+1,4)+offset;
    PExp_xw(:,k) = X(idx,1);
    if PExp_w(k) > 1e-7
        v(:,k) = PExp_xw(:,k) / PExp_w(k);
        p(k) = PExp_w(k); 
    else % if zero or negative
        if PExp_w(k) < -1e-7
            PExp_w
            error('PExp_w should not be negative')
        end
        v(:,k) = zeros(4,1);
        p(k) = 0;
    end
end
%% create list 
p = p / sum(p); % make into a valid distribution
K = N; % number of hypotheses in the list (we consider N hypotheses here)
% list = randsample(N,K,true,p);
list = [1:N];
v_list = v(:,list);
q_list = normalize_cols( v_list );
R_est_list = zeros(3,3,K);
for i=1:K
    qi = q_list(:,i);
    R_est_list(:,:,i) = quat2rot(qi);
end

timeSLIDE = toc(tstart);

%% case of multiple objects in the data
% visualize
q_gt_list = zeros(4,size(R_gt_list,3));
for j = 1:size(R_gt_list,3)
    % get ground truth quaternion for that object
    q_gt_j = rotm2quat(R_gt_list(:,:,j))';
    q_gt_list(:,j)= [q_gt_j(2:4);  q_gt_j(1)]; % put scalar part last
end
q_gt_list
q_list

%% compute minimum errors
R_err_rot_deg_list  = zeros(K,size(R_gt_list,3)); % nr hypotheses vs nr objects 
R_err_norm_list = zeros(K,size(R_gt_list,3)); % nr hypotheses vs nr objects  
for j = 1:size(R_gt_list,3) % nr objects 
    for i=1:K % hypotheses 
        qi = q_list(:,i); % i-th hypothesis
        R_err_rot_deg_list(i,j) = getAngularError(quat2rot(qi),R_gt_list(:,:,j));
        R_err_norm_list(i,j) = norm( vec(quat2rot(qi)) - vec(R_gt_list(:,:,j)) ); 
    end
    fprintf('object nr: %d, R_err_j: %3.2g[deg].\n',j,min(R_err_rot_deg_list(:,j)));
end

%% get error for desired (first) object
R_errs_obj1 = R_err_rot_deg_list(:,1); % error for first object vs hypotheses
[minVal,minInd] = min(R_errs_obj1);
R_err = minVal;
q = q_list(:,minInd);
R_est = quat2rot(q);
fprintf('f_sdp: %3.4e, f_est: %3.4e, R_err: %3.2g[deg].\n',f_sdp,f_est,R_err);
fprintf('Relative suboptimality: %3.2e.\n',subopt);

%% *********************************************
%% *********************************************
%% helper function
function f_est = cost_org(a,b,betasq,q)
R         = quat2rot(q);
residuals = sum((b - R*a).^2,1);
f_est     = sum( min(residuals(:),betasq) );
end

function P = getP()
P = [1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1;
   0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0;
   0, 0, 1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 0, 0;
   0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0;
   -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1;
   0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0;
   0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0;
   0, 0, 0, -1, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, 0;
   -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
end

function [A,b] = maxResidualConstraints(v1,v2,barc2)
P=getP();
P=sparse(P);
N = size(v1,2);
extras = N+1;
n = extras + 4*(N+1);
A = zeros(n,n,N);
b = zeros(N,1);
for k=1:N
    % w_i (||b_i||^2+||a_i||^2 - 2*b_i'*R_i*a_i) \leq barc2
    % w_i (||b_i||^2+||a_i||^2):
    idx = blkIndices(k+1,4) + extras;
    ck = ( v1(:,k)'*v1(:,k)+v2(:,k)'*v2(:,k) );
    A(idx,idx,k) = + ck*eye(4);
    
    % - 2*b_i'*R_i*a_i 
    P_k = reshape(P'*reshape(v2(:,k)*v1(:,k)',[9,1]),[4,4]);
    A([1:4]+ extras,idx,k) = -P_k;
    A(idx,[1:4]+ extras,k) = -P_k;
    
    b(k) = barc2;
end
end

function [Aalpha,balpha] = fixedNumberInliers(N,alpha) % sum w = alpha * N
n = N+1 + 4*(N+1);
Aalpha=zeros(n,n);
Aalpha(1,2:N+1) = ones(1,N);
balpha = alpha * N;
Aalpha = sparse(Aalpha);
end

function B = normalize_cols(A)
mag = sum(A.^2,1).^(0.5);
B   = A./mag;
end

function unitTest(R_gt,theta_gt,N,alpha,v1,v2,betasq,A,b,Aalpha,balpha,isAdversarial)
    %% build X
    q = rotm2quat(R_gt); q = [q(2:4),q(1)]';
    w = (theta_gt + ones(size(theta_gt))) ./ 2;
    v = [[1;w] ; kron([1;w],q)];
    X = v*v';
    extras = N+1;
    q_ind = [1:4]+extras;
    n = extras + 4*(N+1);
    P=getP();
    
    %% check Aalpha,balpha
    if isAdversarial == 0 % in this case theta_gt is not reliable
        assert( abs( trace(Aalpha * X) - balpha )<1e-7 )
    end
    
    %% check A,b
    for k=1:N
        actual = trace( A(:,:,k) * X );
        expected = w(k) * norm(v2(:,k) - R_gt * v1(:,k) )^2; 
        
        assert(abs( actual - expected )<1e-6)
        assert(abs( b(k) - betasq )<1e-6)
    end
        
    %% check constraints
    assert ( abs( trace(X(q_ind,q_ind)) - 1 ) <= 1e-6 ) ;
    assert ( abs( X(1,1) - 1 ) <= 1e-6 );
    for k=1:N
        % top left diagonal (w_i^2 = w_i)
        assert ( norm( X(k+1,k+1) - X(1,k+1) ) <= 1e-6 );
        % bottom right block diagonal (q_i q_i' = w_i^2 q q' = w_i q q' = q q_i)
        idx = blkIndices(k+1,4)+extras;
        assert ( norm( X(idx,idx) - X(q_ind,idx) ) <= 1e-6 );
        % intermediate entries (q_i' = w_i q')
        assert ( norm( X(1,idx) - X(k+1,q_ind) ) <= 1e-6 );
        % intermediate entries (w_i q_i' = w_i^2 q = w_i q = q_i)
        assert ( norm( X(k+1,idx) - X(1,idx) ) <= 1e-6 );
    end
    for k1=1:N
        for k2=k1+1:N+1
            % symmetric off diagonal
            idx1 = blkIndices(k1,4)+extras;
            idx2 = blkIndices(k2,4)+extras;
            assert ( norm( X(idx1,idx2) - X(idx1,idx2)' ) <= 1e-6 );
            % tie theta_i * theta_j to trace of q_i * q_j' (including theta_i vs q * q_i
            assert ( norm( X(k1,k2) - trace( X(idx1,idx2) ) ) <= 1e-6 );
        end
    end
end

