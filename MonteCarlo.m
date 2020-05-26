% Clear memory
clear
close all
clc

rng('Default')
%rng(11); % seedcontrol

%% Monte Carlo / Bootstrap

% Bestimme Anzahl der Pfade
replications     = 1 ; 

% Dementsprechend passe die Spaltenanzahl an
Ncol            = replications;

% Zeithorizont (wird später erhöht)
N       = 50; %30                     % simulation periods excl. transient periods              

% Time horizon
cut             = 1;               % determines the cut of the transient periods
T               = N + (cut-1);      % simulation periods total

%% Parameter Variations
n_spv_vector    = [0.01, 0.1, 0.2];   % Anteil der verbrieften Kredite an Gesamtkrediten (default [0.01, 0.2, 0.5, 0.9])
% For the paper, we run the following scenarios:
% First, we assume a high number of banks (N_B=100) and a relatively high loss absorption of the banking sector (85 %). For this setting, 
% we analyse the effect of different securitization intensities (0.01,0.5,0.9). We then
% change the loss absorption in both directions (75 % and 95 %) and
% evaluate real economic consequences for a securitization intensity of
% 0.1. Finally, we repeat the last analysis for a highly concentrated banking
% sector (N_B=15) and a high loss absorption (95 %).
% n_spv_vector    = 0.05;              % Anteil SPV-Kredit-Verbriefungen
if replications == 1                   % No sensitivity analysis wrt. to the securitization intensity possible 
    n_parameter = 1;                   % Define Benchmark case
else                                   % Sensitivity analysis wrt. to the securitization intensity possible
    n_parameter = numel(n_spv_vector); % Hence, take vector.
end

%% Include numerical values of parameters and exogenous variables
Parameters;                         % Call script Parameters
 
%% Loop that separates the model into regimes by varying securitization intensities
for q = 1: n_parameter              % Loop over securitization intensities
    % Clear some variables
    clearvars ma mi d check ME div N_spvA N_spv div.firms_spvA MC.spv_bankrupt; 
if replications == 1
   n_spv  = n_spv_vector(:,2);      % take the middle entry of the corresponding vector
else
    n_spv  = n_spv_vector(:,q);     % determines the size of the SPV
end
    N_spv  = ceil(n_spv*N_C);       % there are no half firms

    N_spvA   = N_spv;               % determines the size of SPV-A (SPV-B is only implicit)
%% Simulations
% Initialize the variables needed for the monte carlo experiment. The
% folder "MC._" contains the results of the Monte Carlo simulation
MC.g                    = zeros(N,Ncol);        % realized growth rate of the economy (MC1)
MC.gd                   = zeros(N,Ncol);        % expected growth rate of the economy (MC2)
MC.y                    = zeros(N,Ncol);        % realized output (MC3)
MC.y_d                  = zeros(N,Ncol);        % expected output (MC4)
MC.n_bankrupt_i         = zeros(1,Ncol);        % number of insolvent firms (MC5)
MC.n_bankrupt_k         = zeros(1,Ncol);        % number of insolvent banks (MC6)
MC.bankrupt_i           = zeros(N,N_C);         % shows the flag of insolvent firms and the time of insolvency (MC7)          
MC.bankrupt_k           = zeros(N,N_B);         % shows the flag of insolvent banks and the time of insolvency (MC8)
MC.pi_spv               = zeros(N,Ncol);        % profits of the SPV (MC9)
MC.bankrupt_spvA        = zeros(1,N_spvA);      % shows the flag of insolvent firms in the securitized portfolio (MC10)
MC.firms_SPV            = zeros(1,N_spvA);      % shows which corporate loans got securitized (MC11)
MC.SI                   = zeros(N,2);           % checks if ssving is equal to investment (MC12)
MC.CB_solvent           = zeros(N,2);           % checks if assets of solvent firms are equal to liabilities (MC13)
MC.BB_solvent           = zeros(N,2);           % checks if assets of solvent banks are equal to liabilities (MC14)
MC.CB_insolvent         = zeros(N,2);           % checks if assets of insolvent firms are equal to liabilities (MC15)
MC.BB_insolvent         = zeros(N,2);           % checks if assets of insolvent firms are equal to liabilities (MC16)
MC.LR_i                 = zeros(N,N_C);         % leverage ratio of firms (MC17)
MC.ER_k                 = zeros(N,N_B);         % equity ratio of banks (MC18)
MC.bankrupt_firms       = zeros(N,Ncol);        % share of insolvent firms over time (MC19)
MC.bankrupt_banks       = zeros(N,Ncol);        % share of insolvent banks over time (MC20)
MC.spvA                 = zeros(N,Ncol);        % total assets of SPV (MC21)
MC.rmd                  = zeros(N,Ncol);        % !!! contrary to the paper interests on deposits are not constant (MC22)
MC.r_l                  = zeros(N,Ncol);        % lending rate (MC23)
MC.l_i_spv              = zeros(N,N_C);         % shows individual assets of SPV (allows for identifying firm flag) (MC24)
MC.m_d                  = zeros(N,Ncol);        % deposits (MC25)
MC.c                    = zeros(N,Ncol);        % consumption (MC26)
MC.pe_eb                = zeros(N,Ncol);        % equity price of banks (MC27)
MC.pe_ec                = zeros(N,Ncol);        % equity price of firms (MC28)
MC.growth               = zeros(N,Ncol);        % GDP growth rate

% Start simulating the model for "Ncol" replications
    for j = 1:Ncol
%% Load the entire model set-up       
        Model;     
%% Summarize the results and generated time series' in the MC-folder   
        MC.growth(:,j)          = ma.growth(cut:end,:);
        MC.g(:,j)               = ma.g(cut:end,:);                      %MC1
        MC.gd(:,j)              = d.g(cut:end,:);                       %MC2
        MC.y(:,j)               = ma.y(cut:end,:);                      %MC3
        MC.y_d(:,j)             = ma.y_d(cut:end,:);                    %MC4
        MC.n_bankrupt_i(:,j)    = numel(find(mi.bankrupt_i(end,:)==1)); %MC5
        MC.n_bankrupt_k(:,j)    = numel(find(mi.bankrupt_k(end,:)==1)); %MC6
        MC.bankrupt_i(:,:,j)    = mi.bankrupt_i(cut:end,:);             %MC7
        MC.bankrupt_k(:,:,j)    = mi.bankrupt_k(cut:end,:);             %MC8
        MC.pi_spv(:,j)          = ma.pi_spvA(cut:end,:);                %MC9           
        div.bankrupt_spvA_length= length(div.bankrupt_spvA);            %MC10
        MC.bankrupt_spvA(:,1:div.bankrupt_spvA_length,j) = div.bankrupt_spvA;
        MC.firms_SPV(end,:,j)   = div.firms_spvA(end,:);                %MC11
        
%%   Save consistency checks
        MC.SI(:,:,j)            = check.SI(cut:end,:);                  %MC12
        MC.CB_solvent(:,:,j)    = check.CB_solvent(cut:end,:);          %MC13
        MC.BB_solvent(:,:,j)    = check.BB_solvent(cut:end,:);          %MC14
        MC.CB_insolvent(:,:,j)  = check.CB_insolvent(cut:end,:);        %MC15
        MC.BB_insolvent(:,:,j)  = check.BB_insolvent(cut:end,:);        %MC16
        
%%   Continue to summarize results     
        MC.LR_i(:,:,j)          = mi.LR_i(cut:end,:);                   %MC17
        MC.ER_k(:,:,j)          = mi.ER_k(cut:end,:);                   %MC18   
        MC.bankrupt_firms(:,j)  = ma.bankrupt_firms(cut:end,:);         %MC19         
        MC.bankrupt_banks(:,j)  = ma.bankrupt_banks(cut:end,:);         %MC20
        MC.spvA(:,j)            = ma.spvA(cut:end,:);                   %MC21  
        MC.rmd(:,j)             = MC.r_l(:,j) .*MC.m_d(:,j);            %MC22
        MC.r_l(:,j)             = ma.r_l(cut:end,:);                    %MC23
        MC.l_i_spv(:,:,j)       = mi.l_i_spvA(cut:end,:);               %MC24
        MC.m_d(:,j)             = ma.m_d_solvent(cut:end,:);            %MC25
        MC.c(:,j)               = ma.c(cut:end,:);                      %MC26
        MC.pe_eb(:,j)           = ma.pe_eb_solvent(cut:end,:);          %MC27
        MC.pe_ec(:,j)           = ma.pe_ec_solvent(cut:end,:);          %MC28
        
%%   Additional variables not initialized        
        MC.m_dk(:,:,j)          = mi.m_d(cut:end,:);                        %Deposits allocated to single banks                       
        MC.Delta_pe_eb(:,j)     = MC.pe_eb(2:end,j)- MC.pe_eb(1:end-1,j);   %Change of (average) bank equity prices
        MC.Delta_pe_ec(:,j)     = MC.pe_ec(2:end,j)- MC.pe_ec(1:end-1,j);   %Change of (average) firm equity prices


    end
      

 %% Bank Bankruptcies   
        MC.bankrupt_bank_average(:,1)   =  mean(MC.bankrupt_banks,2);       %Average share of bank insolvencies over the replications
        MC.bankrupt_bank_median(:,1)    =  nanmedian(MC.bankrupt_banks,2);  %Median correspondingly 


%% Memory
% Create memories for the different parameter/securitization regimes

ConfidenceBands;                                                            %Does not display outcome at the moment
% "Results does only work with replication >1"
figures; % linked to "Results.m"                                            %Does not display outcome at the moment
% % % % % 
        if q == 1 %Requires 'Results.m'                                     %Save Monte Carlo paths depending on the securitization intensity
            save MC1_workspace1000.mat;
%             print -depsc growthMC1;
        elseif q == 2
            save MC2_workspace1000.mat;
%             print -depsc growthMC2;
        elseif q == 3
            save MC3_workspace1000.mat;
%             print -depsc growthMC3;
        end       
% % % % % save MC_WS.mat stat;



end

%% Anmerkungen:
% Die Monte-Carlo Simulation speichert für jede Parametervariation den
% gesamten 'Workspace' ab, d.h. MC1 beinhaltet die Simulationsergebnisse 
% für den kleinsten Parameterwert von phi_spvA. MC2 den zweitkleinsten, usw...
% Die Datei 'MC_WS' beinhaltet die 'Statistiken' aller Durchläufe.

% In dem Ordner 'MC' des jeweiligen workspaces befinden sich die Werte 
% auserwählter Variablen über jede Schleifeniteration. Diesbezüglich habe
% ich allerdings die ersten 10 Beobachtungen abgeschnitten damit das
% Einschwingungsverhalten nicht Bestandteil unserer Bewertungen sind.

%% Sonstiges
% % figure
% % plot(MC.rmd,'Linewidth',1.5)
% % hold on
% % plot(MC.Delta_pe_eb, 'Linewidth',1.5)
% % % hold on
% % % plot(MC.Delta_pe_ec, 'Linewidth',1.5)
% % hold off
