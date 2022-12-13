%% Test a posteriori bounds in robust rotation search (Wahba Problem)
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
Nset = [50] % number of measurements
outrateSet = [0 0.1 0.2 0.3 0.4 0.5]; % outlier rate
nrTests = 100; % number of runs for each configuration
isAdversarial = 0;

%% allocate matrix to store results
results_true_error = zeros(length(Nset),length(outrateSet),nrTests); 
results_boundN = zeros(length(Nset),length(outrateSet),nrTests); 
results_bound5 = zeros(length(Nset),length(outrateSet),nrTests); 
results_eta = zeros(length(Nset),length(outrateSet),nrTests); 

%% generate random problems and test certifiable hyper-contractivity
for ind_N = 1:length(Nset) % general, but we test for a single N here
    N = Nset(ind_N); % number of measurements
    for ind_outrate = 1:length(outrateSet) 
        outrate = outrateSet(ind_outrate); % outlier rate
        for test = 1:nrTests
            fprintf('=== test %d of %d, with N=%d, and outlier rate = %g ===\n',test,nrTests,N,outrate)
            [problem, R_gt] = createWahbaProblemData(N,outrate,isAdversarial);
            [R_stride,theta_stride,eta] = solveWahba(problem.a,problem.b,problem.betasq,...
                R_gt,... % R_gt only used for performance evaluation
                sdpnalpath,manoptpath,stridepath); % paths to solvers
            
            % compute error and bounds
            true_error = norm( R_gt(:) - R_stride(:) );
            J = []; % set of correctly selected inliers
            for i=1:N
                if theta_stride(i) == problem.theta_gt(i) && theta_stride(i) == 1
                    J = [J i]; 
                end
            end
            
            boundN = computeAposterioriBoundWahba_v2(J, problem.At, sqrt(problem.betasq));
            bound5 = computeAposterioriBoundWahba(theta_stride, problem.At, sqrt(problem.betasq), 5);
            
            % store results
            results_true_error(ind_N,ind_outrate,test) = true_error;
            results_boundN(ind_N,ind_outrate,test) = boundN;
            results_bound5(ind_N,ind_outrate,test) = bound5;
            results_eta(ind_N,ind_outrate,test) = eta;
            
            fprintf('+++ eta=%g, true_error = %g <= %g (boundN) <=? %g (bound5) +++\n',eta,true_error,boundN,bound5)
        end 
    end
end

results_true_error = squeeze(results_true_error);
results_boundN = squeeze(results_boundN);
results_bound5 = squeeze(results_bound5);
results_eta = squeeze(results_eta);

%% plot actual errors vs bounds
figure;
meanErrorTrue = mean(results_true_error,2);
meanBoundN = mean(results_boundN,2);
meanBound5 = mean(results_bound5,2);
semilogy(outrateSet,meanErrorTrue,'-k','linewidth',4);
hold on; grid on
semilogy(outrateSet,meanBoundN,'-r','linewidth',4);
semilogy(outrateSet,meanBound5,'-b','linewidth',4);
semilogy(outrateSet,2*sqrt(3) * ones(size(outrateSet)),'-.k','linewidth',4);
for test = 1:nrTests
  lh = semilogy(outrateSet,results_true_error(:,test),'-k','linewidth',2);
  lh.Color(4)=0.2;
  lh = semilogy(outrateSet,results_boundN(:,test),'-r','linewidth',2); 
  lh.Color(4)=0.2;
  lh = semilogy(outrateSet,results_bound5(:,test),'-b','linewidth',2);
  lh.Color(4)=0.2;
end
semilogy(outrateSet,meanErrorTrue,'-k','linewidth',4);
semilogy(outrateSet,meanBoundN,'-r','linewidth',4);
semilogy(outrateSet,meanBound5,'-b','linewidth',4);
semilogy(outrateSet,2*sqrt(3) * ones(size(outrateSet)),'-.k','linewidth',4);
ylim([1e-4 5e2])
xlabel('Outlier rate ($\beta$)','interpreter','latex')
ylabel('Estimation error ($\|x-x^\circ\|_2$)','interpreter','latex')
legend('actual error','bound-J','bound-5','trivial bound','Location', 'Best')
set(gca,'fontsize', 18)

%% check tightness
figure
results_tightness = results_eta < 1e-7;
tightInstances = 100* sum(results_tightness,2) / nrTests; 
bar(outrateSet,tightInstances)

save log_results_experiment1





