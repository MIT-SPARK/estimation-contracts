%% Visualize bounds in Theorem 11
%% Paper: "Estimation Contracts for Outlier-Robust Geometric Perception"
%% Luca Carlone, Nov 7, 2022

%% plot coefficients C1 and C2:
close all
clc
N = 50;
Ct_pow_t_set = [1 2 4 6];

k = 4;
for j=1:length(Ct_pow_t_set)
    t = k/2;
    Ct_pow_t = Ct_pow_t_set(j); % Ct^t;
    betaMax = ( 1/(Ct_pow_t*2^(3*k-1)) )^(1 / (k/2 - 1)) - 1e-6 % minus some numerical margin
    beta = linspace(0,betaMax,100000);
    C1 = zeros(length(beta),1);
    C2 = zeros(length(beta),1);
    for i=1:length(beta)
        betai = beta(i);
        C1(i) = C1_pow_2_over_k(betai, Ct_pow_t, k);
        C2(i) = C2_pow_2_over_k(betai, Ct_pow_t, k);
    end
    
    %% plot C1
    figure(1)
    semilogy(beta,C1,'-r','linewidth',2)
    hold on
    grid on
    ylimits = [2*10^-3 40];
    xlabel('Outlier rate $\beta$','interpreter','latex')
    ylabel('$C_1(4,\beta)^{1/2}$','interpreter','latex')
    plot([betaMax betaMax],ylimits,'--k','linewidth',1)
    s = sprintf('  C(t)^t = %d  ',Ct_pow_t);
    switch Ct_pow_t
        case 1
            text(beta(end),C1(end)+6,s,'HorizontalAlignment','right')
        case 2
            text(beta(end),C1(end)+4,s,'HorizontalAlignment','left')
        case 4
            text(beta(end),C1(end)+2,s,'HorizontalAlignment','left')
        case 6
            text(beta(end),C1(end)+2,s,'HorizontalAlignment','right')
    end
    ylim(ylimits)
    xlim([-1e-5  5*1e-4])
    
    %% plot C2
    figure(2)
    semilogy(beta,C2,'-b','linewidth',2)
    hold on
    grid on
    ylimits = [1*10^-3 20];
    xlabel('Outlier rate $\beta$','interpreter','latex')
    ylabel('$C_2(4,\beta)^{1/2}$','interpreter','latex')
    plot([betaMax betaMax],ylimits,'--k','linewidth',1)
    ylim([10^-5 10^4])
    s = sprintf(' C(t)^t = %d ',Ct_pow_t);
    switch Ct_pow_t
        case 1
            text(beta(end),C2(end)+1,s,'HorizontalAlignment','right')
        case 2
            text(beta(end),C2(end)+1,s,'HorizontalAlignment','left')
        case 4
            text(beta(end),C2(end)+1,s,'HorizontalAlignment','left')
        case 6
            text(beta(end),C2(end)+1,s,'HorizontalAlignment','right')
    end
    ylim(ylimits)
    xlim([-1e-5  5*1e-4])
end

%% test for increasing k
Ct_pow_t = 6; % Ct^t;
for k=[4 6 8 10]
    t = k/2;
    betaMax = ( 1/(Ct_pow_t*2^(3*k-1)) )^(1 / (k/2 - 1)) - 1e-6 % minus some numerical margin
    beta = linspace(0,betaMax,100000);
    C1 = zeros(length(beta),1);
    C2 = zeros(length(beta),1);
    for i=1:length(beta)
        betai = beta(i);
        C1(i) = C1_pow_2_over_k(betai, Ct_pow_t, k);
        C2(i) = C2_pow_2_over_k(betai, Ct_pow_t, k);
    end
    
    %% plot C1
    figure(3)
    semilogy(beta,C1,'-r','linewidth',2)
    ylimits = [7*1e-5 30];  
    hold on
    grid on
    xlabel('Outlier rate $\beta$','interpreter','latex')
    ylabel('$C_1(k,\beta)^{2/k}$','interpreter','latex')
    plot([betaMax betaMax],ylimits,'--k','linewidth',1)
    s = sprintf('  k = %d  ', k);
    switch k
        case 4
            text(beta(end),C1(end)+5,s,'HorizontalAlignment','left')
        case 6
            text(beta(end),C1(end)+5,s,'HorizontalAlignment','left')
        case 8
            text(beta(end),C1(end)+5,s,'HorizontalAlignment','left')
        case 10
            text(beta(end),C1(end)+5,s,'HorizontalAlignment','right')
    end
    ylim(ylimits)
    xlim([-1e-4  4.5*1e-3])
    
    %% plot C2
    figure(4)
    semilogy(beta,C2,'-b','linewidth',2)
    ylimits = [3*1e-5 10]; 
    hold on
    grid on
    xlabel('Outlier rate $\beta$','interpreter','latex')
    ylabel('$C_2(k,\beta)^{2/k}$','interpreter','latex')
    plot([betaMax betaMax],ylimits,'--k','linewidth',1)
    ylim([10^-5 10^4])
    s = sprintf(' k = %d ',k);
    switch k
        case 4
            text(beta(end),C2(end)+1,s,'HorizontalAlignment','left')
        case 6
            text(beta(end),C2(end)+1,s,'HorizontalAlignment','left')
        case 8
            text(beta(end),C2(end)+1,s,'HorizontalAlignment','left')
        case 10
            text(beta(end),C2(end)+1,s,'HorizontalAlignment','right')
    end
    ylim(ylimits)
    xlim([-1e-4  4.5*1e-3])
end

function value = C1_pow_2_over_k(betai, Ct_pow_t, k)
value = (betai^(k/2-1) * Ct_pow_t * 2^( 3*k-1 ) ) / ...
    ( 1 - betai^(k/2-1) * Ct_pow_t * 2^( 3*k-1 ) );
value = value^(2/k);
end

function value = C2_pow_2_over_k(betai, Ct_pow_t, k)
value = ((2*betai)^(k/2-1) * ( 2^k + Ct_pow_t * 2^( 2*k ) ) ) / ...
            ( 1 - betai^(k/2-1) * Ct_pow_t * 2^( 3*k-1 ) );
value = value^(2/k);
end