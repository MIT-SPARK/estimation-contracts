%% Visualize functions involved in certifiable anti-concentration (definition  in Theorem 11
%% Paper: "Estimation Contracts for Outlier-Robust Geometric Perception"
%% Luca Carlone, Nov 7, 2022

clc
clear
close all
restoredefaultpath
addSpecificPaths

%% set parameters for Monte Carlo runs
Nset = [50] % number of measurements
outrateSet = [0.30]; % outlier rate
nrTests = 5; % number of runs for each configuration
ylimits = [-0.05 1.5]; % for visualization

%% visualize certifiable anti-concentration
figure(1)
ax1 = axes();
for ind_N = 1:length(Nset) % general, but we test for a single N here
    N = Nset(ind_N); % number of measurements
    for ind_outrate = 1:length(outrateSet)
        outrate = outrateSet(ind_outrate); % outlier rate
        [problem, R_gt] = createWahbaProblemData(N,outrate); % only used to get the beta
        for test = 1:nrTests
            fprintf('=== test %d of %d, with N=%d, and outlier rate = %g ===\n',test,nrTests,N,outrate)
            % [problem, R_gt] = createWahbaProblemData(N,outrate);
            % N = problem.N;
            %inliersIds = setdiff([1:N], problem.outliersIds); % anti-concentration only involves the inliers
            %nrInliers = N - problem.nrOutliers;
            alpha = 1-outrate; % inlier rate
            barc = sqrt(problem.betasq); % maximum error for inliers
            
            Mx = sqrt(3); % bound on ||x||
            M = 2*Mx; % bound on || x - x_gt ||
            
            maxSigmaA = 1;
            x = linspace(0,M*maxSigmaA,1000);
            
            %% for bound to be informative
            % alpha * eta / 2 + 2 (1-alpha)/alpha < 2
            eta = (2/alpha) * (2 - 2* (1-alpha)/(alpha))
            delta = 2*barc;
            C = alpha^2 * eta^2 * (1-2*barc)^2 / (32*barc);
            C_delta_Msq = C*delta*M^2;
            
            %% create polynomial
            plotStr = '-k';
            switch test
                case 1
                    c2 = -1;
                    %plotStr = '-r';
                    indText = round(length(x)/4.5);
                    xposText = 0.4;
                    yposText = 0.1;
                case 2
                    c2 = -0.3;
                    indText = round(length(x)/2);
                    %plotStr = '-r';
                    xposText = 1.55;
                    yposText = 0.1;
                case 3
                    c2 = -0.1;
                    indText = round(length(x)/1.2);
                    xposText = 3;
                    yposText = 0.08;
                case 4
                    c2 = -0.05;
                    indText = round(length(x)/1.2);
                    xposText = 2.9;
                    yposText = 0.4;
                case 5
                    c2 = -0.01;
                    indText = round(length(x)/2);
                    %plotStr = '-r';
                    xposText = 2.9;
                    yposText = 0.9;
                otherwise
                    error('exceeded maximum number of cases')
            end
            s = sprintf('  $c_2$ = %g  ', c2);
            
            y = (1 + c2 * x.^2).^2;
            
            %% create upper bound
            for i=1:length(x)
                yup(i) = C_delta_Msq / x(i)^2;
            end
            
            %% create lower bounds
            for i=1:length(x)
                if x(i) < delta
                    yd1(i) = (1-delta)^2;
                else
                    yd1(i) = 0;
                end
            end
            for i=1:length(x)
                if x(i) < delta
                    yd2(i) = (1+delta)^2;
                else
                    yd2(i) = 2*ylimits(2);
                end
            end
            
            %% plot
            plot(x,yup,'-.r','linewidth',2);
            hold on; grid on
            plot(x,y,plotStr,'linewidth',2);
            plot(x,yd1,':k','linewidth',2);
            %plot(x,yd2,':k','linewidth',2);
            text(xposText,yposText,s,'HorizontalAlignment','left', 'Interpreter', 'latex','FontSize',16)
            ylim(ylimits)
            
            xlabel('$\| A_i  v\|$', 'Interpreter', 'latex')
            ylabel('$p^2(\|A_i v\|)$ and constraints', 'Interpreter', 'latex')
            legend('$C \delta M^2 / \|v\|^2$','$p^2(\|A_i v\|)$','$(1 - \delta)^2$','Location', 'ne','interpreter','latex')
            set(gca,'fontsize', 18)
        end
    end
end
text(1.7,1.6,'$\|v\|$','HorizontalAlignment','left','color','r', 'Interpreter', 'latex','FontSize',18)
ax2 = axes('Position', get(ax1,'Position'), ...
    'XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k');
ax2.YAxis.Visible = 'off';
ax2.XLim = ax1.XLim;
ax2.XTick = ax1.XTick;
ax2.XColor = 'r';


%% plot increasingly stricter thresholds
% set parameters for Monte Carlo runs
Nset = [50] % number of measurements
outrateSet = [0.2 0.3 0.4 0.45]; % outlier rate
nrTests = 1; % number of runs for each configuration
ylimits = [-0.05 1.5]; % for visualization

figure(2)
ax1 = axes();
% generate random problems and visualize certifiable anti-contraction
for ind_N = 1:length(Nset) % general, but we test for a single N here
    N = Nset(ind_N); % number of measurements
    for ind_outrate = 1:length(outrateSet)
        outrate = outrateSet(ind_outrate); % outlier rate
        [problem, R_gt] = createWahbaProblemData(N,outrate); % only used to get the beta
        
        alpha = 1-outrate; % inlier rate
        barc = sqrt(problem.betasq); % maximum error for inliers
        
        Mx = sqrt(3); % bound on ||x||
        M = 2*Mx; % bound on || x - x_gt ||
        
        maxSigmaA = 1;
        x = linspace(0,M*maxSigmaA,1000);
        
        %% for bound to be informative
        % alpha * eta / 2 + 2 (1-alpha)/alpha < 2
        alpha
        eta = (2/alpha) * (2 - 2* (1-alpha)/(alpha))
        delta = 2*barc;
        C = alpha^2 * eta^2 * (1-2*barc)^2 / (32*barc);
        C_delta_Msq = C*delta*M^2;
        
        %% create polynomial
        plotStr = '-k';
        switch ind_outrate
            case 1
                xposText = 3;
                yposText = 0.66;
            case 2
                xposText = 3;
                yposText = 0.39;
            case 3
                xposText = 3;
                yposText = 0.15;
            case 4
                xposText = 1.5;
                yposText = 0.18;
            otherwise
                error('exceeded maximum number of cases')
        end
        
        c2 = -0.1;
        indText = round(2/3*length(x));
        s = sprintf('  $\\eta$ = %.2f  ', eta);
        
        y = (1 + c2 * x.^2).^2;
        
        %% create upper bound
        for i=1:length(x)
            yup(i) = C_delta_Msq / x(i)^2;
        end
        
        %% create lower bounds
        for i=1:length(x)
            if x(i) < delta
                yd1(i) = 1-delta;
            else
                yd1(i) = 0;
            end
        end
        for i=1:length(x)
            if x(i) < delta
                yd2(i) = 1+delta;
            else
                yd2(i) = 2*ylimits(2);
            end
        end
        
        
        %% plot
        plot(x,yup,'-.r','linewidth',2);
        hold on; grid on
        plot(x,y,plotStr,'linewidth',2);
        plot(x,yd1,':k','linewidth',2);
        %plot(x,yd2,':k','linewidth',2);
        text(xposText,yposText,s,'HorizontalAlignment','left', 'Interpreter', 'latex','FontSize',16)
        ylim(ylimits)
        
        xlabel('$\| A_i  v\|$', 'Interpreter', 'latex')
        ylabel('$p^2(\|A_i v\|)$ and constraints', 'Interpreter', 'latex')
        legend('$C \delta M^2 / \|v\|^2$','$p^2(\|A_i v\|)$','$(1 - \delta)^2$','Location', 'ne','interpreter','latex')
        set(gca,'fontsize', 18)
        
        etaSet(ind_outrate) = eta;
    end
end
text(1.7,1.6,'$\|v\|$','HorizontalAlignment','left','color','r', 'Interpreter', 'latex','FontSize',18)
ax2 = axes('Position', get(ax1,'Position'), ...
    'XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k');
ax2.YAxis.Visible = 'off';
ax2.XLim = ax1.XLim;
ax2.XTick = ax1.XTick;
ax2.XColor = 'r';

%% plot bounds
% set parameters for Monte Carlo runs
N = [50] % number of measurements
nrTests = 4; % number of runs for each configuration
ylimits = [-0.05 1.5]; % for visualization

alpha = [0.5:0.01:1];

outrate = outrateSet(ind_outrate); % outlier rate
[problem, R_gt] = createWahbaProblemData(N,outrate); % only used to get the beta
for etaInd = 1:length(etaSet)
    eta = etaSet(etaInd);
    fprintf('=== test %d of %d, with N=%d, and outlier rate = %g ===\n',test,nrTests,N,outrate)
    for i = 1:length(alpha)
        ai = alpha(i);
        b(i) = Mx * (ai * eta) / 2 + 2 * Mx * (1-ai) / ai;
    end
    
    switch etaInd
        case 1
            xposText = 0.93;
            yposText = 3.1;
        case 2
            xposText = 0.93;
            yposText = 2.7;
        case 3
            xposText = 0.93;
            yposText = 2.15;
        case 4
            xposText = 0.93;
            yposText = 1.4;
        otherwise
            error('exceeded maximum number of cases')
    end
    indText = round(4/5*length(alpha));
    s = sprintf('  $\\eta$ = %.2f  ', eta)
    
    %% plot
    figure(3)
    semilogy(alpha,b,'-b','linewidth',2)
    hold on
    semilogy([alpha],2 * Mx * ones(size(alpha)),'--r','linewidth',2)
    text(xposText,yposText,s,'HorizontalAlignment','left', 'Interpreter', 'latex','fontsize',16)
    grid on
    box off
    % xticks([0 0.2 0.4 0.5 0.6 0.8 1])
    ylabel('Estimation error bounds','interpreter','latex')
    xlabel('Inlier rate $\alpha$','interpreter','latex')
    legend('proposed', 'trivial bound','interpreter','latex')
    ax = gca;
    ax.FontSize = 18;
end

%% visualize certifiable anti-concentration for higher degree polynomials
Nset = [50] % number of measurements
outrateSet = [0.45]; % outlier rate
nrTests = 3; % number of runs for each configuration
ylimits = [-0.05 1.5]; % for visualization
degreeSet = [2 4 8];

figure(4)
ax1 = axes();
for ind_N = 1:length(Nset) % general, but we test for a single N here
    N = Nset(ind_N); % number of measurements
    for ind_outrate = 1:length(outrateSet)
        outrate = outrateSet(ind_outrate); % outlier rate
        [problem, R_gt] = createWahbaProblemData(N,outrate); % only used to get the beta
        for test = 1:nrTests
            degree = degreeSet(test)
            fprintf('=== test %d of %d, with N=%d, and outlier rate = %g ===\n',test,nrTests,N,outrate)
            alpha = 1-outrate; % inlier rate
            barc = sqrt(problem.betasq); % maximum error for inliers
            
            Mx = sqrt(3); % bound on ||x||
            M = 2*Mx; % bound on || x - x_gt ||
            
            maxSigmaA = 1;
            x = linspace(0,M*maxSigmaA,1000);
            
            %% for bound to be informative
            eta = (2/alpha) * (2 - 2* (1-alpha)/(alpha))
            delta = 2*barc;
            C = alpha^2 * eta^2 * (1-2*barc)^2 / (32*barc);
            C_delta_Msq = C*delta*M^2;
            
            %% fit coefficients:
            xplot= [0:0.001:1.2*M];
            yindicator = (abs(xplot)<=delta);
            nrPoints = length(xplot);
            
            %% degree
            
            A = zeros(nrPoints,1);
            b = yindicator(:)-ones(nrPoints,1);
            for i=1:nrPoints
                for d = 2:2:degree
                    A(i,d/2) = [xplot(i)^(d)];
                end
            end
            c = inv(A' * A) * A' * b

            %% create polynomial
            plotStr = '-k';
            switch test
                case 1
                    y = (1 + c(1) * x.^2).^2;
                    xposText = 1.7;
                    yposText = 0.55;
                    
                case 2
                    y = (1 + c(1) * x.^2 + c(2) * x.^4).^2;
                    indText = round(length(x)/2);
                    %plotStr = '-r';
                    xposText = 1.25;
                    yposText = 0.4;
                case 3
                    y = (1 + c(1) * x.^2 + c(2) * x.^4 + c(3) * x.^6 + c(4) * x.^8).^2;
                    xposText = 0.38;
                    yposText = 0.2;
                otherwise
                    error('exceeded maximum number of cases')
            end
            s = sprintf('  deg = %d  ', degree);
            
            %% create upper bound
            for i=1:length(x)
                yup(i) = C_delta_Msq / x(i)^2;
            end
            
            %% create lower bounds
            for i=1:length(x)
                if x(i) < delta
                    yd1(i) = (1-delta)^2;
                else
                    yd1(i) = 0;
                end
            end
            for i=1:length(x)
                if x(i) < delta
                    yd2(i) = (1+delta)^2;
                else
                    yd2(i) = 2*ylimits(2);
                end
            end
            
            %% plot
            plot(x,yup,'-.r','linewidth',2);
            hold on; grid on
            plot(x,y,plotStr,'linewidth',2);
            plot(x,yd1,':k','linewidth',2);
            %plot(x,yd2,':k','linewidth',2);
            text(xposText,yposText,s,'HorizontalAlignment','left', 'Interpreter', 'latex','FontSize',16)
            ylim(ylimits)
            
            xlabel('$\| A_i  v\|$', 'Interpreter', 'latex')
            ylabel('$p^2(\|A_i v\|)$ and constraints', 'Interpreter', 'latex')
            legend('$C \delta M^2 / \|v\|^2$','$p^2(\|A_i v\|)$','$(1 - \delta)^2$','Location', 'ne','interpreter','latex')
            set(gca,'fontsize', 18)
        end
    end
end
text(1.7,1.6,'$\|v\|$','HorizontalAlignment','left','color','r', 'Interpreter', 'latex','FontSize',18)
ax2 = axes('Position', get(ax1,'Position'), ...
    'XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k');
ax2.YAxis.Visible = 'off';
ax2.XLim = ax1.XLim;
ax2.XTick = ax1.XTick;
ax2.XColor = 'r';
