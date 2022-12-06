function [isHyperContractive,Q,Z,timeTest] = testHypercontractivity(At,Ct_pow_t_input)

%% create polynomial in definition of hyper-contractivity
N = size(At,3); % nr measurements

pvar v1 v2 v3 v4 v5 v6 v7 v8 v9
v = [v1;v2;v3;v4;v5;v6;v7;v8;v9]; 
for i=1:N
    normsq_Ai_v(i) = sum( ( At(:,:,i) * v ).^2 ); % || A_i' v ||^2
end

t = 2;
% Ct = N^( (t-1) / t )

if nargin<2
    Ct_pow_t = N^( (t-1) ) % justified in the paper
else
    warning('using user-specified Ct')
    Ct_pow_t = Ct_pow_t_input;
end

%% polynomial in eq (26) of the paper
p_hyper = Ct_pow_t * ( (1/N) * sum( normsq_Ai_v ) )^t ...
    - (1/N) * sum( normsq_Ai_v.^t ); % degree 4 polynomial for t=2

%% check if it is sos using SOSTOOLS
tstart = tic;
[Q,Z] = findsos(p_hyper); 
timeTest = toc(tstart);

%% if sos, the decomposition Q, Z is not empty
isHyperContractive = ~isempty(Q) && ~isempty(Z); % true if p_hyper is sos