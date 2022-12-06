function [R_stride,theta_stride,eta] = solveWahba(a,b,betasq,R_gt,sdpnalpath,manoptpath,stridepath)
%% solve SDP relaxation using STRIDE solver, scalable and faster
%% see Yang et al. An Inexact Projected Gradient Method with Rounding and Lifting 
%% by Nonlinear Programming for Solving Rank-One Semidefinite Relaxation of Polynomial Optimization
%% https://arxiv.org/abs/2105.14033
%% Solve using STRIDE

% generate standard SDP data (this is different than CVX)
SDP = QUASAR_Problem(a,b,betasq);

% set parameters for STRIDE
options.pgdStepSize     = 10; % step size, default 10
options.maxiterPGD      = 5; % maximum outer iterations for STRIDE, default 5-10
options.SDPNALpath      = sdpnalpath; % provide path to SDPNAL
options.tolADMM         = 5e-6; % tolerance for warmstart, decrease this parameter for a better warmstart (but takes more time)
options.tolPGD          = 1e-8; % tolerance on KKT residual of the SDP
options.lbfgseps        = false;
% provide implementation to the local search method
options.rrOpt           = 1:3; % round the leading 3 eigenvectors to generate hypotheses
options.rrFunName       = 'local_search_quasar'; % name of the .m file that implements the local search

% Primal initialization
[R_gnc,info_gnc]    = GNC_Wahba(a,b,betasq,1.4);
q_gnc               = rotm2quat(R_gnc); q_gnc = [q_gnc(2:4),q_gnc(1)]';
v_gnc               = kron([1;info_gnc.theta_gnc],q_gnc);
X0                  = {v_gnc*v_gnc'};

% call STRIDE
addpath(genpath(manoptpath));
addpath(genpath(sdpnalpath));
addpath(genpath(stridepath));
[outPGD,Xopt,yopt,Sopt] = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, X0, options);

if size(R_gt,3) > 1
    R_gt = R_gt(:,:,1);
end 

infostride              = get_performance_quasar(Xopt,yopt,Sopt,SDP,R_gt);
R_stride = infostride.R; 
theta_stride = infostride.theta(:);
eta = infostride.Rs; % suboptimality gap

%% *********************************************
%% *********************************************
%% helper function
function Q_cost = cost(v1,v2,barc2)
P=[1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1;
   0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0;
   0, 0, 1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 0, 0;
   0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0;
   -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1;
   0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0;
   0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0;
   0, 0, 0, -1, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, 0;
   -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
P=sparse(P);
N = size(v1,2);
Npm = 4*N + 4;
Q_1=zeros(Npm,Npm);
for k=1:N
    idx = 4+blkIndices(k,4);
    P_k = reshape(P'*reshape(v2(:,k)*v1(:,k)',[9,1]),[4,4]);
    ck = 0.5 * ( v1(:,k)'*v1(:,k)+v2(:,k)'*v2(:,k) - barc2 );
    Q_1((1:4),idx) = Q_1((1:4),idx)-0.5*P_k+ck/2*eye(4);
    Q_1(idx,(1:4)) = Q_1(idx,(1:4))-0.5*P_k+ck/2*eye(4);
end

Q_2=zeros(Npm,Npm);
for k=1:N
%     idx = 4+blkIndices(k,4);
    idx = blkIndices(1,4);
    P_k = reshape(P'*reshape(v2(:,k)*v1(:,k)',[9,1]),[4,4]);
    ck = 0.5 * ( v1(:,k)'*v1(:,k)+v2(:,k)'*v2(:,k) + barc2 );
    Q_2(idx,idx) = Q_2(idx,idx) - P_k + ck*eye(4);
end
Q_cost=Q_1+Q_2;
Q_cost=sparse(Q_cost);

Q_cost=Q_cost/barc2;
end

function B = normalize_cols(A)
mag = sum(A.^2,1).^(0.5);
B   = A./mag;
end

function f_est = cost_org(a,b,betasq,q)
R         = quat2rot(q);
residuals = sum((b - R*a).^2,1)/betasq;
f_est     = sum( min(residuals(:),1) );
end

end
