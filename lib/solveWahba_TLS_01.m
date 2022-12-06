%% solve SDP TLS relaxation using cvx (with w in {0;1})

v1=problem.a;
v2=problem.b;
betasq=problem.betasq;

N = size(v1,2); % nr vectors (measurements)
% X = [q\tran, q_1\tran, q_2\tran, ..., q_N\tran] has size n = 4+4*N
n = 4*N + 4;

%% solve SDP relaxation
C = cost(v1,v2,betasq);
unitTest(R_gt,problem.theta_gt(:),N,v1,v2,betasq,C);
q_ind = [1:4];

cvx_begin % sdp
% cvx_solver mosek
variable X(n,n) symmetric
minimize( trace(C * X) ) %
subject to
X == semidefinite(n);
trace(X(q_ind,q_ind)) == 1;
for k=1:N
    % w_i^2 = w_i
    idx = blkIndices(k+1,4);
    X(idx,idx) == X(q_ind,idx);
end
for k1=1:N
    for k2=k1+1:N+1
        % symmetric off diagonal
        idx1 = blkIndices(k1,4);
        idx2 = blkIndices(k2,4);
        X(idx1,idx2) == X(idx1,idx2)';
    end
end
cvx_end

% extract solution
f_sdp   = cvx_optval; % lower bound
eig(X)

%% ROUNDING
[V,~]   = sorteig(X);
v       = V(:,1);
q       = normalize_cols( v(1:4) );
w_beforeRounding   = zeros(N,1);
for i = 1:N
    w_beforeRounding(i) = round(q'*v(blkIndices(i+1,4)));
end
w = round(w_beforeRounding);

R_err   = getAngularError(quat2rot(q),R_gt);
f_est   = cost_org(v1,v2,betasq,q); % upper bound
subopt  = abs(f_est - f_sdp) / (1+abs(f_est)+abs(f_sdp));

eta = subopt;
R_est = quat2rot(q);
fprintf('f_sdp: %3.4e, f_est: %3.4e, R_err: %3.2e[deg].\n',f_sdp,f_est,R_err);
fprintf('Relative suboptimality: %3.2e.\n',subopt);

%% *********************************************
%% *********************************************
%% helper function
function Q_cost = cost(v1,v2,barc2)
P=getP();
P=sparse(P);
N = size(v1,2);
n = 4*N + 4;
Q_cost=zeros(n,n);
for k=1:N
    % w_i (||b_i||^2+||a_i||^2 - barcq) + barcq - 2*b_i'*R_i*a_i
    % barcq:
    Q_cost([1:4],[1:4]) = Q_cost([1:4],[1:4]) + barc2*eye(4);
    
    % w_i (||b_i||^2+||a_i||^2 - barcq):
    idx = blkIndices(k+1,4);
    ck = ( v1(:,k)'*v1(:,k)+v2(:,k)'*v2(:,k) - barc2 );
    Q_cost(idx,idx) = Q_cost(idx,idx) + ck*eye(4);
    
    % - 2*b_i'*R_i*a_i 
    P_k = reshape(P'*reshape(v2(:,k)*v1(:,k)',[9,1]),[4,4]);
    Q_cost([1:4],idx) = Q_cost([1:4],idx)-P_k;
    Q_cost(idx,[1:4]) = Q_cost(idx,[1:4])-P_k;
end
Q_cost=sparse(Q_cost);
end

function B = normalize_cols(A)
mag = sum(A.^2,1).^(0.5);
B   = A./mag;
end

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

function unitTest(R_gt,theta_gt,N,v1,v2,barc2,C)
    %% build X
    q = rotm2quat(R_gt); q = [q(2:4),q(1)]';
    w = (theta_gt + ones(size(theta_gt))) ./ 2;
    v = kron([1;w],q);
    X = v*v';
    n = 4*(N+1);
    P=getP();
    
    %% check A,b
    for k=1:N
        Q_cost = zeros(n,n);
        % w_i (||b_i||^2+||a_i||^2 - barcq) + barcq - 2*b_i'*R_i*a_i
        % barcq:
        Q_cost([1:4],[1:4]) = Q_cost([1:4],[1:4]) + barc2*eye(4);
        
        % w_i (||b_i||^2+||a_i||^2 - barcq):
        idx = blkIndices(k+1,4);
        ck = ( v1(:,k)'*v1(:,k)+v2(:,k)'*v2(:,k) - barc2 );
        Q_cost(idx,idx) = Q_cost(idx,idx) + ck*eye(4);
        
        % - 2*b_i'*R_i*a_i
        P_k = reshape(P'*reshape(v2(:,k)*v1(:,k)',[9,1]),[4,4]);
        Q_cost([1:4],idx) = Q_cost([1:4],idx)-P_k;
        Q_cost(idx,[1:4]) = Q_cost(idx,[1:4])-P_k;
        
        cost_expected = w(k)*norm(v2(:,k) - R_gt*v1(:,k))^2 + (1-w(k))*barc2;
        
        assert(abs( trace(Q_cost * X) - cost_expected )<1e-6)
    end
    
    assert(abs( trace(C * X) - cost_org(v1,v2,barc2,q) )<1e-6)
end

