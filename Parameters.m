%% General Numerical Values 
N_C             = 300;                 % Number of firms
N_B             = 50; %30,50,70        % Number of banks
par.searchcost  = 5;                   % Searchcost (scaling-multiplicator)
ER_T            = 0.05; %0.09          % Thresholdt equity ratio
LR_T            = 0.9;                 % Threshold leverage ratio
 
% Date of securitization
t_spv           = 20; %20

%% Concerning initialisation
par.theta_spvA  = 0.95;   %0.95         % Scaling parameter SPVA-return           
par.theta_c     = 0.20;   %0.34         % Determines initial endowment of firm's equities
par.theta_b     = 0.07;  %0.115         % Determines initial endowment of bank's equities
par.theta_r     = 0.85;   %0.9 %0.85    % weight of the mark-up in the loan rate

%% Market Entry
par.eta      =  5; %5%8%6               % max. number of recapitalised bank's per period
par.zeta     =  2;  %12%1%3%6           % Duration of the insolvency proceedings
par.kappa    = 0.2; %0.2%0.5%0.05       % size of recapitalisation

%% Firms
par.xi_c     = 0.2;                    % fraction of distributed profits
par.delta    = 0.000;                  % capital depreciation
par.gamma    = 0.5; %0.7,0.85,1.0 %0.5 % percentage of spv credit that is swapped for deposits instead of bank securities
par.gamma_1  = 0.03; %0.02;            % coefficient of the rate of capacity utilization
par.gamma_2  = 0.025;                  % coefficient of the profit share
par.gamma_3  = 0.3;                    % coefficient of the profit rate

par.lambda   = 0.01;                   % Adaptive Expectation Parameter 
par.mu_c     = 0;                      % Mean consumption expectation shock
par.sigma_c  = 1.0;                    % Std. consumption expectation shock
%% Households
par.omega_0  = 0.6;                    % wage share %0.6 default
par.alpha_1  = 0.6;  %0.6% 0.12        % propensity to consume out of dispodable income 
%par.alpha_2  = 0.02; % 0.1             % autonomous consumption
par.alpha_3  = 0.055;%0.055;            % propensity to consume out of deposits 
par.alpha_4  = 0.045;%0.045;            % propensity to consume out of corporate equities and tranche A
par.alpha_5  = 0.03; %0.03;%0.01;       % propensity to consume out of bank equities

%% Banking Sector
par.r_l0     = 0.04;                   % long-term average loan rate
par.r_md0    = par.theta_r*par.r_l0;   % rate on deposits

par.rho      = 0.06;                   % Scaling Parameter linking the leverage ratio to the loan rate
par.rho_ER   = 0.010;                  % Scaling Parameter linking the bank's equity ratio to the loan rate
LR_bar       = 0.50;                   % Bank's benchmark leverage ratio

% Properties of the stochasti terms in the laon rate
par.sigma_rl = 0.007;                  % std. shock in loan rate's mark-up
par.mean_rl  = 0;                      % mean of the normal distr. shock

%% Initialisierung
K_0         = 200;                   % Anfangsausstattung Kapitalstock
y_0         = 0.5*K_0;               % Anfangsausstattung BIP
c_0         = 0.7*y_0;               % Anfangsausstattung Konsum
p_ei0       = 1;                     % Aktienpreise nicht-finanzielle Unternehmen (eigentlich Initialisierung)
p_ek0       = 1;                     % Aktienpreise Banken (eigentlich Initialisierung)
M           = 5;                     % Anzahl potenzieller Geschäftspartner

%%
par.rho_mu  = 0.0;                 % Auslastungsschock Ar Koeffizient
par.sigma_mu = 0.1;
%% Firms asset prices
par.sig0 = 0.05;
par.sig1 = 0.05;
par.sig2 = 0.05;





save BaselineCalibration.mat;