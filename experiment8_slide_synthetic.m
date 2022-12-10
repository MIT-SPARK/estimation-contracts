%% Test list decodable regression in robust rotation search (Wahba Problem)
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
outrateSet = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]; % outlier rate
nrTests = 10; % number of runs for each configuration
% isAdversarial = 0; recoverAllHypotheses = 0; % Fig. 9
%isAdversarial = 1; recoverAllHypotheses = 0; % Fig. 10
 isAdversarial = 1; recoverAllHypotheses = 1; % Fig. 11
doPlotRotations = 0;

%% parameters to visualize rotations out of SLIDE (Fig. 12)
% Nset = [50] % number of measurements
% outrateSet = 0.8 % [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]; % outlier rate
% nrTests = 1; % number of runs for each configuration
% isAdversarial = 3;
% doPlotRotations = 1;

%% allocate matrix to store results
results_TLS_R_err_ang_deg_obj1 = zeros(length(Nset),length(outrateSet),nrTests); 
results_TLS_R_err_diff_norm_obj1 = zeros(length(Nset),length(outrateSet),nrTests); 
results_TLS_eta = zeros(length(Nset),length(outrateSet),nrTests); 
results_TLS_R_err_ang_deg_obj2 = zeros(length(Nset),length(outrateSet),nrTests); 
results_TLS_R_err_diff_norm_obj2 = zeros(length(Nset),length(outrateSet),nrTests); 
results_TLS_fsdp = zeros(length(Nset),length(outrateSet),nrTests); 

results_SLIDE_R_err_ang_deg_obj1 = zeros(length(Nset),length(outrateSet),nrTests); 
results_SLIDE_R_err_diff_norm_obj1 = zeros(length(Nset),length(outrateSet),nrTests); 
results_SLIDE_eta_norm = zeros(length(Nset),length(outrateSet),nrTests); 
results_SLIDE_eta = zeros(length(Nset),length(outrateSet),nrTests); 
results_SLIDE_R_err_ang_deg_obj2 = zeros(length(Nset),length(outrateSet),nrTests);
results_SLIDE_R_err_diff_norm_obj2 = zeros(length(Nset),length(outrateSet),nrTests); 
results_SLIDE_fsdp = zeros(length(Nset),length(outrateSet),nrTests); 
results_SLIDE_time = zeros(length(Nset),length(outrateSet),nrTests); 

%% generate random problems and test certifiable hyper-contractivity
for ind_N = 1:length(Nset) % general, but we test for a single N here
    N = Nset(ind_N); % number of measurements
    for ind_outrate = 1:length(outrateSet) 
        outrate = outrateSet(ind_outrate); % outlier rate
        for test = 1:nrTests
            fprintf('=== test %d of %d, with N=%d, and outlier rate = %g ===\n',test,nrTests,N,outrate)
            [problem, R_gt, R_gt_list] = createWahbaProblemData(N,outrate,isAdversarial);
            
            if doPlotRotations==1
                for i=1:size(R_gt_list,3)
                    my_ref_frame(R_gt_list(:,:,i),zeros(3,1),['r','g','b'],0.75,0.05,0.5)
                end
            end

            %% Different TLS formulation, for testing
            solveWahba_TLS_m1p1
            
            %% store results single hypotheses
            results_TLS_fsdp(ind_N,ind_outrate,test) = f_sdp;
            results_TLS_R_err_ang_deg_obj1(ind_N,ind_outrate,test) = R_err;
            R_err_diff_norm = norm( R_gt(:) - R_est(:) );
            results_TLS_R_err_diff_norm_obj1(ind_N,ind_outrate,test) = R_err_diff_norm;
            results_TLS_eta(ind_N,ind_outrate,test) = eta;
            fprintf('+++ eta=%g +++\n',eta)
            fprintf('+++ TLS obj 1: rot err =%.3g [deg], norm error=%.3g +++\n',R_err,R_err_diff_norm)
            if isAdversarial > 0 % compute error for object 2
                for j=2:size(R_gt_list,3)
                    R_err   = getAngularError(quat2rot(q),R_gt_list(:,:,j));
                    results_TLS_R_err_ang_deg_obj2(ind_N,ind_outrate,test) = R_err;
                    R_err_diff_norm = norm( vec(quat2rot(q)) - vec(R_gt_list(:,:,j)) );
                    if j == 2 % we only log the second one, but display all
                        results_TLS_R_err_diff_norm_obj2(ind_N,ind_outrate,test) = R_err_diff_norm;
                    end
                    fprintf('+++ TLS obj %d: rot err =%.3g [deg], norm error=%.3g +++\n',j, R_err,R_err_diff_norm)
                end
            end
            clear eta R_err q R_est R_err_diff_norm
            fprintf('========================================\n')
            
            %% List decodable
            solveWahba_listDecodable_01
            
            %% store results list decodable
            results_SLIDE_time(ind_N,ind_outrate,test) = timeSLIDE;
            results_SLIDE_fsdp(ind_N,ind_outrate,test) = f_sdp;
            results_SLIDE_R_err_ang_deg_obj1(ind_N,ind_outrate,test) = R_err;
            R_err_diff_norm = norm( R_gt(:) - R_est(:) );
            results_SLIDE_R_err_diff_norm_obj1(ind_N,ind_outrate,test) = R_err_diff_norm;
            if R_err_diff_norm > 2*sqrt(3)
                R_err_diff_norm
                error('error exceeds trivial upper bound?')
            end
            results_SLIDE_eta_norm(ind_N,ind_outrate,test) = eta_norm;
            results_SLIDE_eta(ind_N,ind_outrate,test) = eta;
            fprintf('+++ eta=%g +++\n',eta)
            fprintf('+++ SLIDE obj 1: rot err =%.3g [deg], norm error=%.3g +++\n',R_err,R_err_diff_norm)
            if isAdversarial > 0
                assert(isnan(R_err) || abs( R_err - min(R_err_rot_deg_list(:,1)) )<1e-7 )
                results_SLIDE_R_err_ang_deg_obj2(ind_N,ind_outrate,test) = min(R_err_rot_deg_list(:,2));
                results_SLIDE_R_err_diff_norm_obj2(ind_N,ind_outrate,test) = min(R_err_norm_list(:,2));
                fprintf('+++ SLIDE obj 2: rot err =%.3g [deg], norm error=%.3g +++\n',min(R_err_rot_deg_list(:,2)),min(R_err_norm_list(:,2)))
            end  
            
            if doPlotRotations==1
               plotRotations(R_gt_list, R_est_list) 
            end
            clear eta R_err q R_est R_err_diff_norm R_err_list v_list q_list R_est_list X
        end
    end
end

save log_results_experiment8_isAdversarial_1_recoverAll_1
 
%% plots
% clear all 
% close all
% clc
% load log_results_experiment7_isAdversarial_true
s = [1:length(outrateSet)]; 
outrateSet = outrateSet(s);
margin = 0.02;

results_TLS_R_err_ang_deg_obj1 = squeeze(results_TLS_R_err_ang_deg_obj1);
results_TLS_R_err_ang_deg_obj1 = results_TLS_R_err_ang_deg_obj1(s,:);
results_TLS_R_err_diff_norm_obj1 = squeeze(results_TLS_R_err_diff_norm_obj1);
results_TLS_R_err_diff_norm_obj1 = results_TLS_R_err_diff_norm_obj1(s,:);
results_TLS_eta = squeeze(results_TLS_eta);
results_TLS_eta = results_TLS_eta(s,:);
results_TLS_R_err_ang_deg_obj2 = squeeze(results_TLS_R_err_ang_deg_obj2);
results_TLS_R_err_ang_deg_obj2 = results_TLS_R_err_ang_deg_obj2(s,:);
results_TLS_R_err_diff_norm_obj2 = squeeze(results_TLS_R_err_diff_norm_obj2);
results_TLS_R_err_diff_norm_obj2 = results_TLS_R_err_diff_norm_obj2(s,:);

results_SLIDE_R_err_ang_deg_obj1 = squeeze(results_SLIDE_R_err_ang_deg_obj1); 
results_SLIDE_R_err_ang_deg_obj1 = results_SLIDE_R_err_ang_deg_obj1(s,:);
results_SLIDE_R_err_diff_norm_obj1 = squeeze(results_SLIDE_R_err_diff_norm_obj1);
results_SLIDE_R_err_diff_norm_obj1 = results_SLIDE_R_err_diff_norm_obj1(s,:);
results_SLIDE_eta = squeeze(results_SLIDE_eta);
results_SLIDE_eta = results_SLIDE_eta(s,:);
results_SLIDE_R_err_ang_deg_obj2 = squeeze(results_SLIDE_R_err_ang_deg_obj2);
results_SLIDE_R_err_ang_deg_obj2 = results_SLIDE_R_err_ang_deg_obj2(s,:);
results_SLIDE_R_err_diff_norm_obj2 = squeeze(results_SLIDE_R_err_diff_norm_obj2);
results_SLIDE_R_err_diff_norm_obj2 = results_SLIDE_R_err_diff_norm_obj2(s,:);
results_SLIDE_time = squeeze(results_SLIDE_time);

% plot errors (norm of difference)
figure
plot(outrateSet, mean(results_TLS_R_err_diff_norm_obj1,2) ,'-b','linewidth',4);
hold on; grid on
plot(outrateSet, mean(results_SLIDE_R_err_diff_norm_obj1,2) ,'-r','linewidth',4);
plot_distribution_prctile(outrateSet,results_TLS_R_err_diff_norm_obj1','color',[0 0 1],'Prctile',[25 50 75 90]);
plot_distribution_prctile(outrateSet,results_SLIDE_R_err_diff_norm_obj1','color',[1 0 0],'Prctile',[25 50 75 90]);
plot(outrateSet, mean(results_TLS_R_err_diff_norm_obj1,2) ,'-b','linewidth',4);
plot(outrateSet, mean(results_SLIDE_R_err_diff_norm_obj1,2) ,'-r','linewidth',4);
% plot(outrateSet,2*sqrt(3) * ones(size(outrateSet)),'-.k','linewidth',4);
xlabel('Outlier rate ($\beta$)','interpreter','latex')
ylabel('Estimation error ($\|x-x^\circ\|_2$)','interpreter','latex')
% legend('TLS','SLIDE','trivial bound','Location', 'Best')
legend('TLS','SLIDE','Location', 'northwest')
set(gca,'fontsize', 18)
xlim([outrateSet(1)-margin outrateSet(end)+margin])
yl = ylim; ylim([-0.05*(yl(2)-yl(1)) yl(2)])

% plot errors (rotation errors)
figure
plot(outrateSet, mean(results_TLS_R_err_ang_deg_obj1,2) ,'-b','linewidth',4);
hold on; grid on
plot(outrateSet, mean(results_SLIDE_R_err_ang_deg_obj1,2) ,'-r','linewidth',4);
plot_distribution_prctile(outrateSet,results_TLS_R_err_ang_deg_obj1','color',[0 0 1],'Prctile',[25 50 75 90]);
plot_distribution_prctile(outrateSet,results_SLIDE_R_err_ang_deg_obj1','color',[1 0 0],'Prctile',[25 50 75 90]);
plot(outrateSet, mean(results_TLS_R_err_ang_deg_obj1,2) ,'-b','linewidth',4);
plot(outrateSet, mean(results_SLIDE_R_err_ang_deg_obj1,2) ,'-r','linewidth',4);
xlabel('Outlier rate ($\beta$)','interpreter','latex')
ylabel('Rotation error [deg]','interpreter','latex')
legend('TLS','SLIDE','Location', 'northwest')
set(gca,'fontsize', 18)
xlim([outrateSet(1)-margin outrateSet(end)+margin])
yl = ylim; ylim([-0.05*(yl(2)-yl(1)) yl(2)])

% relaxation gap
figure
semilogy(outrateSet, mean(results_TLS_eta,2) ,'-b','linewidth',4);
hold on;
semilogy(outrateSet, mean(results_SLIDE_eta,2) ,'-r','linewidth',4);
plot_distribution_prctile(outrateSet,results_TLS_eta','color',[0 0 1],'Prctile',[25 50 75 90]);
plot_distribution_prctile(outrateSet,results_SLIDE_eta','color',[1 0 0],'Prctile',[25 50 75 90]);
semilogy(outrateSet, mean(results_TLS_eta,2) ,'-b','linewidth',4);
semilogy(outrateSet, mean(results_SLIDE_eta,2) ,'-r','linewidth',4);
xlabel('Outlier rate ($\beta$)','interpreter','latex')
ylabel('Relaxation gap','interpreter','latex')
legend('TLS','SLIDE','Location', 'Best')
set(gca,'fontsize', 18)
grid on
xlim([outrateSet(1)-margin outrateSet(end)+margin])

if isAdversarial
    if exist('recoverAllHypotheses','var') ==0 || recoverAllHypotheses == 0
        % plot errors (rotation errors) for each object
        figure
        plot(outrateSet, mean(results_TLS_R_err_ang_deg_obj1,2) ,'-b','linewidth',4);
        hold on; grid on
        plot(outrateSet, mean(results_TLS_R_err_ang_deg_obj2,2) ,':b','linewidth',4);
        plot(outrateSet, mean(results_SLIDE_R_err_ang_deg_obj1,2) ,'-r','linewidth',4);
        plot(outrateSet, mean(results_SLIDE_R_err_ang_deg_obj2,2) ,':r','linewidth',4);
        plot_distribution_prctile(outrateSet(5:end),results_SLIDE_R_err_ang_deg_obj2(5:end,:)','color',[1 0 0],'Prctile',[25 50 75 90]);
        plot_distribution_prctile(outrateSet(1:5),results_SLIDE_R_err_ang_deg_obj2(1:5,:)','color',[0.5 0.5 0.5],'Prctile',[25 50 75 90]);
        plot(outrateSet, mean(results_TLS_R_err_ang_deg_obj1,2) ,'-b','linewidth',4);
        plot(outrateSet, mean(results_TLS_R_err_ang_deg_obj2,2) ,':b','linewidth',4);
        plot(outrateSet, mean(results_SLIDE_R_err_ang_deg_obj1,2) ,'-r','linewidth',4);
        plot(outrateSet, mean(results_SLIDE_R_err_ang_deg_obj2,2) ,':r','linewidth',4);
        xlabel('Outlier rate ($\beta$)','interpreter','latex')
        ylabel('Rotation error [deg]','interpreter','latex')
        legend('TLS - rotation 1','TLS - rotation 2','SLIDE - rotation 1','SLIDE - rotation 2','Location', 'best')
        set(gca,'fontsize', 18)
        xticks([0.1:0.1:0.9])
        xlim([outrateSet(1)-margin outrateSet(end)+margin])
        yl = ylim; ylim([-0.05*(yl(2)-yl(1)) yl(2)])
    else
        % plot errors (rotation errors) for each object
        figure
        plot(outrateSet, mean(results_TLS_R_err_ang_deg_obj1,2) ,'-b','linewidth',4);
        hold on; grid on
        plot(outrateSet, mean(results_TLS_R_err_ang_deg_obj2,2) ,':b','linewidth',4);
        plot(outrateSet, mean(results_SLIDE_R_err_ang_deg_obj1,2) ,'-r','linewidth',4);
        plot(outrateSet, mean(results_SLIDE_R_err_ang_deg_obj2,2) ,':r','linewidth',4);
        plot_distribution_prctile(outrateSet(5:end),results_SLIDE_R_err_ang_deg_obj2(5:end,:)','color',[1 0 0],'Prctile',[25 50 75 90]);
        plot_distribution_prctile(outrateSet(1:5),results_SLIDE_R_err_ang_deg_obj2(1:5,:)','color',[1 0 0],'Prctile',[25 50 75 90]);
        %plot_distribution_prctile(outrateSet,results_SLIDE_R_err_ang_deg_obj1','color',[1 0 0],'Prctile',[25 50 75 90]);
        plot(outrateSet, mean(results_TLS_R_err_ang_deg_obj1,2) ,'-b','linewidth',4);
        plot(outrateSet, mean(results_TLS_R_err_ang_deg_obj2,2) ,':b','linewidth',4);
        plot(outrateSet, mean(results_SLIDE_R_err_ang_deg_obj1,2) ,'-r','linewidth',4);
        plot(outrateSet, mean(results_SLIDE_R_err_ang_deg_obj2,2) ,':r','linewidth',4);
        xlabel('Outlier rate ($\beta$)','interpreter','latex')
        ylabel('Rotation error [deg]','interpreter','latex')
        legend('TLS - rotation 1','TLS - rotation 2','SLIDE - rotation 1','SLIDE - rotation 2','Location', 'best')
        set(gca,'fontsize', 18)
        xticks([0.1:0.1:0.9])
        xlim([outrateSet(1)-margin outrateSet(end)+margin])
        yl = ylim; ylim([-0.05*(yl(2)-yl(1)) yl(2)])
    end
end
