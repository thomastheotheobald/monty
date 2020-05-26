
% -------------------------------------------------------------------------
% Initialisierung
% -------------------------------------------------------------------------
initval.omega = par.omega_0;                % wage share (assumed to be constant)
initval.K   = K_0;                          % aggregated capital stock
initval.c   = c_0;                          % aggregated consumption
initval.y   = y_0;                          % output
initval.u   = y_0/K_0;                      % utilization
initval.wb  = initval.omega*initval.y;      % wage bill
initval.e_c = par.theta_c*K_0;              % equity stock firms
initval.p_ec= p_ei0;                        % equity price firms
initval.l_c0= K_0 - initval.e_c*initval.p_ec;    % credit stock (macro)          
initval.p_eb= p_ek0;                        % equity stock banks 

%% Initialize randomly firm's leverage ratio
    mi.LR_i         = zeros(T,N_C);           % leverage ratio of single firms  
    for z= 1:N_C
        mi.LR_i(1,z)= unifrnd(0.5,0.85);       % draws uniformly between the values
    end
    mi.l_i          = zeros(T,N_C);           % credit of single firms
    mi.l_i(1,:)     = mi.LR_i(1,:)*(K_0/N_C); % initialize credit of single firms

%%    
    mi.l_rest                   = zeros(T,N_C);                            % ???
    div.bestbank_i              = zeros(T,N_C);                            % records the firm-bank contracts
    div.N_BB                    = randi(N_B,N_C,1);                        % in t=1, realtionships are random                               
    div.bestbank_i(1,:)         = div.N_BB;                                % returns an uniformly distributed random integer
    ma.l_c                      = zeros(T,1);                              % credit (macro)
    ma.l_c(1)                   = sum(mi.l_i(1,:));                        % initial value
    
    ma.l_cs                     = zeros(T,1);                              % for the ev of the model not necessary !!!
    ma.l_cs(1)                  = ma.l_c(1);                               % "
    initval.l_c0 = ma.l_cs(1);                                             % "
    
    ma.l_c_insolvent            = zeros(T,1);                               % credit of insolvent firms                              
    ma.l_c_solvent              = zeros(T,1);                               % credit of solvent firms
    ma.l_c_solvent(1,1)         = ma.l_cs(1);                               % initial value
    %% Endogenize firm's equities
    mi.e_i                      = zeros(T,N_C);
    mi.e_i(1,:)                 = ((K_0/N_C)-mi.l_i(1,:))/p_ei0;            % firm equity micro level
    mi.var_pei                  = zeros(T,N_C);
         mi.var_pei(1,:) = 1;  
    mi.sigma_pei                = zeros(T,N_C);
    
    %%
    mi.p_ei                     = zeros(T,N_C);                             % equity price firms micro level
    mi.p_ei(1,:)                = 1;          % initial value 1
%     mi.p_ei(1,:)                = ((K_0/N_C)-mi.l_i(1,:))/initval.e_c;          % initial value
        

    ma.p_ec                     = zeros(T,1);                               % equity price firms macro level
    ma.p_ec(1)                  = (sum(mi.p_ei(1,:))*initval.e_c)/(N_C*initval.e_c);  % initial value
    d.p_ec                      = zeros(T,1);                               % production planning equity price
    d.p_ec(1)                   = ma.p_ec(1);                               % initial value
initval.e_b                     = par.theta_b*ma.l_cs(1);                   % initial share stocks of total bank balance
initval.e_k                     = initval.e_b/N_B;                          % stock equity initialisation 
    mi.e_k                      = zeros(T,N_B);                             % bank equities micro level
    mi.e_k(1,:)                 = initval.e_k;                              % initial value
    ma.e_b                      = zeros(T,1);                               % bank equities macro level
    ma.e_b(1)                   = initval.e_b;                              % initial value
    ma.l_bank_insolvent         = zeros(T,1);                               % lending by insolvent banks macro level
    ma.l_bank_solvent           = zeros(T,1);                               % lending by solvent banks macro level
    ma.l_bank_solvent(1,:)      = ma.l_c_solvent(1,1);                      % initial value
    mi.l_k_insolvent            = zeros(T,N_B);                             % lending by insolvent banks micro level
initval.m_d                     = initval.l_c0 - initval.e_b*initval.p_eb;  % deposits (macro)
initval.pi_T                    = y_0 - initval.wb;                         % total profits (macro)
initval.pi_D                    = par.xi_c*initval.pi_T;                    % distributed profits
initval.pi_tilde                = initval.pi_T - initval.pi_D;              % retained earnings
initval.V_h                     = initval.m_d + p_ei0*initval.e_c + p_ek0*initval.e_b;   % Household's Wealth initial value
initval.i                       = par.delta*K_0;                                         % Investment initial value
initval.g                       = par.delta;                                             % growth rate initial value
initval.y_d                     = initval.c;                                             % Disposable Income inital value
ma.spvA                         = zeros(T,1);                               % Securitized loans
div.l_ki_spvA_NB                = zeros(T,N_B);                             % Securitized loans micro level                                      % Kredite der B-Tranche
ma.r_spvA                       = zeros(T,1);                               % Interest rate of securitized loans macro level
div.r_l_ki_NB                   = zeros(T,N_B);                             % Interest rate of securitized loans mirco level
  if N_spvA > 0   
   div.firms_spvA = zeros(T,N_spvA);                                        % Individual firm loans involved in the securitization
  else
   div.firms_spvA = zeros(T,1);                                             % Individual firm loans involved in the securitization
  end
mi.l_i_spvA                     = zeros(T,N_C);                             % Securitized loans micro level 
mi.r_li_spvA                    = zeros(T,N_C);                             % Interest rates micro level
mi.l_i_exspvA                   = zeros(T,N_C);                             % Loans micro level without securitized loans
mi.l_i_exspvA(1,:)              = mi.l_i(1,:);                              % Initial value
ma.pi_spvA                      = zeros(T,1);                               % Profits of special purpose vehicle
div.bankrupt_spvA               = 0;                                        % Flag for insolvent firms involved in securitization
%% ------------------------------------------------------------------------
% Mikro
% -------------------------------------------------------------------------
initval.p_ei        = p_ei0;                                            % initial firm equity price
initval.p_ek        = p_ek0;                                            % initial bank equity price
initval.l_i         = initval.l_c0/N_C;                                 % initial credit demand micro level
initval.Phi         = 1/N_C;                                            % initial size of corporates
initval.pi_i        = initval.Phi*(initval.pi_T)- par.r_l0*initval.l_i; % profits of individual non-financial firms (micro)
initval.pi_tilde_i  = (1-par.xi_c)*initval.pi_i;                        % retained earning of individual firm (micro)
initval.LR_i        = initval.l_i/(initval.l_i + initval.p_ei*initval.e_c);  % Firm's Leverage Ratio
initval.g_i         = initval.g*initval.Phi;
% Banking Sector
%% Variablen Initialisierung
% Makroebene
ma.eq_ch            = zeros(T,1);
ma.K                = zeros(T,1);                                           % capital stock maro level       
ma.K(1)             = K_0;                                                  % initial value
d.K                 = zeros(T,1);                                           % desired capital stock macro level
d.K(1)              = K_0;                                                  % initial value
% -------------------------------------------------------------------------
ma.K_insolvent      = zeros(T,1);                                           % capital stock insolvent firms macro level 
ma.K_solvent        = zeros(T,1);                                           % capital stock solvent firms macro level
ma.K_solvent(1,1)   = K_0;                                                  % initial values
ma.pi_T_solvent     = zeros(T,1);                                           % profits solvent firms macro level
ma.pi_T_solvent(1,1)= initval.pi_T;                                         % initial values
ma.e_c_insolvent    = zeros(T,1);                                           % insolvent firms' equities macro level 
ma.e_c_solvent      = zeros(T,1);                                           % solvent firms' equities macro level
ma.e_c_solvent(1,1) = initval.e_c;                                          % initial Value
ma.pe_ec_insolvent  = zeros(T,1);                                           % insolvent firms' equity price macro level
ma.pe_ec_solvent    = zeros(T,1);                                           % solvent firms' equity price macro level
ma.pe_ec_solvent    = sum(mi.p_ei(1,:)*initval.e_c/N_C);                             % initial value
mi.p_ek             = zeros(T,N_B);                                         % bank eqity prices micro level
mi.p_ek(1,:)        = p_ek0;                                                % initial value
ma.pe_eb_solvent    = zeros(T,1);                       % cross product of equities and equity prices solvent banks
ma.pe_eb_solvent    = sum(mi.p_ek(1,:)*initval.e_k);    % initial value 
ma.e_b_insolvent    = zeros(T,1);                       % insolvent banks' equity macro level
ma.e_b_solvent      = zeros(T,1);                       % solvent banks' equity macro level
ma.e_b_solvent(1,1) = initval.e_b;                      % inital value
ma.pe_eb_insolvent  = zeros(T,1);                       % cross products of equities and equity prices insolvent banks
ma.m_d_insolvent    = zeros(T,1);                       % deposits insolvent banks    
% -------------------------------------------------------------------------
% desired variables stem from production planning
ma.bankrupt_firms   = zeros(T,1);                       % share of insolvent firms macro level
ma.bankrupt_banks   = zeros(T,1);                       % share of insolvent banks macro level
ma.u                = zeros(T,1);                       % capacity utilization (macro)
ma.u(1)             = initval.u;                        % initialize
d.u                 = zeros(T,1);                       % desired capacity utilization (macro)
d.u(1)              = initval.u;                        % initialize
ma.g                = zeros(T,1);                       % growth rate (macro)
ma.g(1)             = initval.g;                        % initialize
d.g                 = zeros(T,1);                       % desired growth rate (macro)
d.g(1)              = initval.g;                        % initialize
ma.y                = zeros(T,1);                       % GDP (macro)
ma.y(1)             = initval.y;                        % initialize
d.y                 = zeros(T,1);                       % desired GDP
d.y(1)              = initval.y;                        % initialize
ma.y_d              = zeros(T,1);                       % household disposable income (macro) 
ma.y_d(1)           = initval.y_d;                      % initialize
ma.c                = zeros(T,1);                       % consumption (macro)
ma.c(1)             = initval.c;                        % initialize
d.c                 = zeros(T,1);                       % desired consumption (macro)
d.c(1)              = initval.c;                        % initialize
ma.i                = zeros(T,1);                       % investment (macro)
ma.i(1)             = initval.i;                        % initialize
d.i                 = zeros(T,1);                       % desired investment (macro)
d.i(1)              = initval.i;                        % initialize
ma.pi_T             = zeros(T,1);                       % total profits (macro)
ma.pi_T(1)          = initval.pi_T;                     % initialize
d.pi_T              = zeros(T,1);                       % desired total profits (macro)
d.pi_T(1)           = initval.pi_T;                     % initialize
ma.pi_D             = zeros(T,1);                       % distributed profits (macro)
ma.pi_D(1)          = initval.pi_D;                     % initialize
d.pi_D              = zeros(T,1);                       % desired distributed profits (macro)
d.pi_D(1)           = initval.pi_D;                     % initialize
ma.pi_tilde         = zeros(T,1);                       % retained earnings (macro)
ma.pi_tilde(1)      = initval.pi_tilde;                 % initialize
d.pi_tilde          = zeros(T,1);                       % desired retained earnings (macro)
d.pi_tilde(1)       = initval.pi_tilde;                 % initialize
ma.V_h              = zeros(T,1);                       % household wealth (macro)
ma.V_h(1)           = initval.V_h;                      % initialize
ma.p_eb             = zeros(T,1);                       % bank equity prices (macro)
ma.p_eb(1)          = p_ek0;                            % initialize
ma.e_c              = zeros(T,1);                       % corporate equities (macro)
ma.e_c(1)           = initval.e_c;                      % initialize
ma.r_l                      = zeros(T,1);               % lending rate (macro)
ma.m_ds                     = zeros(T,1);               % deposit rate (macro)
ma.m_ds(1,:)                = initval.m_d;              % initialize
ma.p_eb_insolvent           = zeros(T,1);               % equity prices of insolvent banks
ma.omega                    = zeros(T,1)+0.7;           % wage share
ma.cons_share               = zeros(T,1);               % consumption share of GDP
ma.cons_share(1,:)          = ma.c(1,:)/ma.y(1,:);      % initialize
ma.prof_rate                = zeros(T,1);               % profit rate
ma.prof_rate(1,:)           = ma.pi_T(1,:)/ma.K(1,:);   % initialize
ma.growth                   = zeros(T,1);               % GDP growth (g = growth of capital stock)
ma.y_d_share                = zeros(T,1);               % share of household disposable income
ma.y_d_share                = ma.y_d(1,:)/ma.y(1,:);    % initialize
% mi. für Mikroebene
mi.Phi_i                    = zeros(T,1);               % individual proportion of the desired variables, similar to market share                                        
mi.Phi_i(1,1)               = 1/N_C;
mi.r_li                     = zeros(T,N_C);             % lending rates (micro)
mi.r_li(1,:)                = par.r_l0;                 % initialize
mi.i_i                      = zeros(T,N_C);             % investment (micro)
mi.i_i(1,:)                 = (par.delta*K_0)/N_C;      % initialize
mi.pi_i                     = zeros(T,N_C);             % total profits (micro)
mi.pi_i(1,:)                = initval.pi_i;             % initialize
mi.pi_iD                    = zeros(T,N_C);             % distributed profits (micro)
mi.pi_iD(1,:)               = par.xi_c*initval.pi_i;    % initialize
mi.pi_tilde_i               = zeros(T,N_C);             % retained earnings (micro)
mi.pi_tilde_i(1,:)          = initval.pi_tilde_i;       % initialize
mi.r_lbar_i                 = zeros(T,1);               % reference lending rate (micro)
mi.r_lbar_i(1,:)            = par.r_l0;                 % initialize
mi.bankrupt_i               = zeros(T,N_C);             % flag for insolvent firms    
%%  Banking Sector   
%  random distribution of the laons (random network)
div.M_k                     = randsample(1:N_C,N_C);    % draw randomly banks for each of the firms (oversized)
initval.digit_k             = floor(N_C/N_B);           % firm-bank ratio (integer) 
div.l                       = mi.l_i(1,div.M_k);        % firms' credit allocation at the beginning of the simulation
mi.l_s_ik                   = zeros(T,N_B);             % credit granted by banks (micro level)
mi.l_s_ik(1,1)              = sum(div.l(1,1:initval.digit_k)); % initialize (first initval.digit_k loans go to first bank)
mi.l2_s_ik                  = zeros(T,N_B);             % data graveyard
mi.md2                      = zeros(T,N_B);             % bank deposits (micro level)
% compute the entire loan volume of digit_k randomly picked firms 
for z =2:N_B
  mi.l_s_ik(1,z)            = sum(div.l(1,(z*initval.digit_k-(initval.digit_k -1)):z*initval.digit_k));
end
% example for z =2: initially bank 2 grants firm loans 5:8, note that this
% fully works when the number of firms is an integer mulitple of the number
% of banks
    if sum(mi.l_i(1,:))             > sum(mi.l_s_ik(1,:))
        mi.l_s_ik(1,randperm(N_B,1))= mi.l_s_ik(1,randperm(N_B,1))...
                                    + sum(mi.l_i(1,:)) - sum(mi.l_s_ik(1,:));
    end
% if this is not the case, allocate the rest
% Bank's Equity Ratio
mi.ER_k(1,:)                = initval.p_ek*initval.e_k./(mi.l_s_ik(1,:));   % bank equity ratios (micro level)  
mi.rest_k                   = zeros(T,N_B);                                 % banks' rest of SPV's balance (should be 0)
mi.r_mk                     = zeros(T,N_B);                                 % deposit rates (micro level)
mi.r_mk(1,:)                = par.r_md0;                                    % initialize
initval.Phi_k               = ma.l_c(1)/N_B;                                % initial credit share of each bank
mi.m_d                      = zeros(T,N_B);                                 % deposits (micro level)
mi.m_d(1,:)                 =  mi.l_s_ik(1,:)-initval.p_ek*initval.e_k;     % initialize so that banks assets are equal to liabilities
mi.pi_k                     = zeros(T,N_B);                                 % total bank profits (micro level)                                        
mi.pi_k(1,:)                = initval.Phi*initval.pi_T;                     % initialize
mi.bankrupt_k               = zeros(T,N_B);                                 % flag for insolvent banks
mi.K_i                      = zeros(T,N_C);                                 % firm capital stock (micro level)
mi.K_i(1,:)                 = K_0/N_C;                                      % initialize
mi.K_i_insolvent            = zeros(T,N_C);                                 % capital stock of insolvent firms (micro level)
mi.transaction_cost         = zeros(T,N_C);                                 % flag times number of insolvent bank drawings
mi.l_i_delta                = zeros(T,N_C);                                 % corporate credit growth (micro level)

%% Market Entries (bank recapitalization)
div.ER_sort                 = zeros(T,N_B);                                 % data graveyard (sort banks by equity ratio)
div.ER_order                = zeros(T,N_B);                                 % data graveyard (bank flag by equity ratio)
ME.count                    = zeros(T,N_B);                                 % counts the number of periods during bank insolvency
ME.bestER                   = zeros(T,N_B);                                 % displays the (absolute) proportion of best banks' equity ratio used for ME
mi.ER2                      = zeros(T,N_B);                                 % datagraveyard 
d.mu                        = zeros(T,1);
