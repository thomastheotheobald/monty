

%% Lade die Anfangsausstattungen (Initialwerte)
Initialisation;


%% Simulations    
 for t =2:T
% disp(t)
%% ------------------------------------------------------------------------
% Macro Level
% -------------------------------------------------------------------------
% Ar Process utilization shock
    d.mu(t,:)   = par.rho_mu*d.mu(t-1,:) + normrnd(0,par.sigma_mu);
% Bank's and Firm's equities    
    ma.e_c(t,:) = ma.e_c_solvent(t-1,:);    % are equal to those of the previous period
    ma.e_b(t,:) = ma.e_b_solvent(t-1,:);    % are equal to those of the previous period

% desired capital stock growth rate (*or before interaction)
    d.g(t,:) =   par.gamma_1*d.u(t-1,:) + ...  % depending on capacity utilization, profit share and profit rate   
        par.gamma_2*(1-ma.omega(t-1,:)) + par.gamma_3*(ma.pi_T(t-1,:)/ma.K(t-1,:));
% desired investment (*or just planned investment)
    d.i(t,:)   = (d.g(t,:)+par.delta)*ma.K(t-1,:); % is equal to multiplicative growth of previous period's capital stock 
% planned stock of capital
    d.K(t,:)   = d.i(t,:) + (1-par.delta)*ma.K(t-1,:); % accumulation of capital by investment
    
% desired stock of firm´ equities
    d.p_ec(t,:) = (1+d.g(t,:))*ma.p_ec(t-1,:);         % equity prices grow in line with economic growth
    d.p_eb(t,:) = (1+d.g(t,:))*ma.p_eb(t-1,:);         % equity prices grow in line with economic growth

% desired consumption 
%     d.c(t,:)   = par.alpha_1*ma.y_d(t-1,:)...             % depends on household disposable income and wealth    
%         +par.alpha_2*ma.V_h(t-1,:);
    d.c(t,:)   = ma.c(t-1,:) + par.lambda*(d.c(t-1,:) - ma.c(t-1,:))...    % Adaptive expectations   
        + normrnd(par.mu_c, par.sigma_c);
% Output and capacity utilization which is associated with the the planned investment and consumption
    d.y(t,:) = d.c(t,:) + d.i(t,:);                     % desired output consist of consumption and investment
    d.u(t,:) = d.y(t,:) / (d.K(t,:)*(1+d.mu(t,:)));     % utilization is defined as a percentage of capital stock
    
% profits
    d.pi_T(t,:)        = (1-ma.omega(t,:))*d.y(t,:);    % total profits are determined by profit share
    d.pi_D(t,:)        = par.xi_c*d.pi_T(t,:);          % distributed profits are determined by payout ratio

         if div.bankrupt_spvA == 0                                              % if spvA exists                           
            d.pi_tilde(t,:) = d.pi_T(t,:) - d.pi_D(t,:)...                      % desired retained earnings are equal to total profits minus distributed
                - sum(mi.r_li(t-1,:).*(mi.l_i(t-1,:)-mi.l_i_spvA(t-1,:))) ...   % minus cost of debt without spv
                - sum(mi.r_li_spvA(t-1,:).*mi.l_i_spvA(t-1,:));                 % minus cost of securitized loans 
         else                                                                   % if spvA is bankrupt
            d.pi_tilde(t,:) = d.pi_T(t,:) - d.pi_D(t,:)...                      % desired retained earnings are equal to total profits minus distributed
                - sum(mi.r_li(t-1,:).*(mi.l_i(t-1,:)-mi.l_i_spvA(t-1,:))) ...   % minus capital cost of debt less the securitized loans
                - sum(mi.r_li_spvA(t-1,:).*mi.l_i_spvA(t-1,:));...               % minus capital cost of securitized loans
%               -sum(mi.r_li_spvA(t-1,div.bankrupt_spvA).*mi.l_i_spvA(t-1,div.bankrupt_spvA));
%               Pausibility ? Programm is running when commented out.
         end
% firms credit demand (desired or before interaction)
    ma.l_c(t,:) = ma.l_c_solvent(t-1,:) + d.i(t,:) - d.pi_tilde(t,:);           % new credit finances rest of investment which is not covered by internal financing    

    
%% Bank-Firm Interaction at the Micro-Level
% Consists of:
% - partner selection mechanism
% - firm sector micro level
% - banking sector at the micro with and without securitization
% - SPV
% - market entry/exit: incl. bank bailout 
Interaction;
%% Macro Level after interaction (Bottom-Up)
       
% Constraint Macro Variables  
% Determine Macro variables by taking into account firms realized
% investment and credit demand. Therefore, planned investment differs from
% actual or realized investment due to restrictions at the micro level.
            ma.K(t,:) = sum(mi.K_i(t,:));                                               % aggregate individual capital stocks
            ma.g(t,:)   = (sum(mi.K_i(t,:)) - sum(mi.K_i(t-1,:)))/sum(mi.K_i(t-1,:));   % determine realized growth rate from capital stock
            ma.i(t,:)   = sum(mi.i_i(t,:));                                             % aggregate individual investment
            ma.K_insolvent(t,:) = sum(mi.K_i(t, mi.bankrupt_i(t,:) ==1));               % aggregate capital stock of insolvent firms
            ma.K_solvent(t,:) = sum(mi.K_i(t, mi.bankrupt_i(t,:) ==0));                 % aggregate capital stock of solvent firms
%%     
% Depositen
            ma.m_d_insolvent(t,:) = sum(mi.m_d(t, mi.bankrupt_k(t,:)==1 ));             % aggregate deposits of insolvent banks          
            ma.m_d_solvent(t,:) =sum(mi.m_d(t, mi.bankrupt_k(t,:)==0 ));                % aggregate deposits of solvent banks
            ma.m_d(t,:) =sum(mi.m_d(t,:));                                              % aggregate deposits of all
% Equities
            ma.p_ec(t,:) = sum(mi.p_ei(t, mi.bankrupt_i(t,:)==0).*mi.e_i(t, mi.bankrupt_i(t,:)==0))/sum(mi.e_i(t, mi.bankrupt_i(t,:)==0));%(1+ma.g(t,:))*ma.p_ec(t-1,:);                                % firm equity prices grow with aggregated growth rate
            ma.e_c_insolvent(t,:) = sum(mi.e_i(t,mi.bankrupt_i(t,:) ==1));              % find equities of insolvent firms        
            ma.e_c_solvent(t,:) = sum(mi.e_i(t,mi.bankrupt_i(t,:) == 0));               %ma.e_c(t-1,:);                                        % equities of solvent firms do not change ?
            ma.pe_ec_insolvent(t,:) = sum(mi.p_ei(t, mi.bankrupt_i(t,:)==1)...
                .*mi.e_i(t,mi.bankrupt_i(t,:) ==1));                          % compute product of equities and equity prices for insolvent firms                      
            ma.pe_ec_solvent(t,:) = sum(mi.p_ei(t, mi.bankrupt_i(t,:)==0).*mi.e_i(t, mi.bankrupt_i(t,:)==0));     % compute product of equities and equity prices for solvent firms        
            ma.e_b_solvent(t,:) = sum(mi.e_k(t,mi.bankrupt_k(t,:) == 0));               % find equity of solvent banks
            ma.e_b_insolvent(t,:) = sum(mi.e_k(t,mi.bankrupt_k(t,:) ==1));              % find equity of insolvent banks
            ma.pe_eb_solvent(t,:) = sum(mi.p_ek(t,mi.bankrupt_k(t,:) == 0).*mi.e_k(t,mi.bankrupt_k(t,:) == 0)); % compute product of equities and equity prices for solvent banks
            ma.pe_eb_insolvent(t,:) = sum(mi.p_ek(t, mi.bankrupt_k(t,:)==1)...          % compute product of equities and equity prices for insolvent banks
                .*mi.e_k(t, mi.bankrupt_k(t,:)==1));  % cross-product          
            ma.p_eb(t,:) = (ma.l_bank_solvent(t,:)-ma.m_d_solvent(t,:))./ma.e_b_solvent(t,:); % compute equity prices of solvent banks in a way that closes the balances
            ma.p_eb_insolvent(t,:) = ma.p_eb_insolvent(t,:)+ ...
                ~isnan((ma.l_bank_insolvent(t,:)-ma.m_d_insolvent(t,:))./ma.e_b_insolvent(t,:)); % compute equity prices of insolvent banks in way that closes the balance

%% Household Sector
%  realized consumption
            ma.c(t,:)   = par.alpha_1*ma.y_d(t-1,:)...                      % consumption depends on disposable income and wealth    
                        + par.alpha_3*ma.m_d(t-1,:)...
                        + par.alpha_4*(ma.pe_ec_solvent(t-1,:)+ma.pe_ec_insolvent(t-1,:)+ma.spvA(t-1,:))...
                        + par.alpha_5*(ma.pe_eb_solvent(t-1,:)+ma.pe_eb_insolvent(t-1,:));
% Consumption depending on todays income               
            ma.y(t,:)   = ma.c(t,:) + ma.i(t,:);                            % gdp consist of consumption and investment
            ma.u(t,:)   = ma.y(t,:) / ma.K(t,:);                            % utilization as a percentage of capital stock
% corporate profits
            ma.pi_T(t,:)= (1-ma.omega(t,:))*ma.y(t,:);                      % profits depend on exogeneously set profit share
            ma.pi_T_solvent(t,:) = sum(mi.pi_i(t,mi.bankrupt_i(t,:)==0));   % profits of solvent firms are aggregated from the micro level: how is this connected to the previous line ?
            ma.pi_D(t,:)= par.xi_c*ma.pi_T(t,:);                            % distributed profits depend on exogeneously set payout ratio
         
         if div.bankrupt_spvA == 0                                              % if the spv exists,
            ma.pi_tilde(t,:) = ma.pi_T(t,:) - ma.pi_D(t,:)...                   % retained earnings are equal to total profits minus distributed profits
                - sum(mi.r_li(t-1,:).*(mi.l_i(t-1,:)-mi.l_i_spvA(t-1,:))) ...   % minus capital cost of debt less the securitized loans
                - sum(mi.r_li_spvA(t-1,:).*mi.l_i_spvA(t-1,:));                 % minus capital cost of the securitized loans
         else
            ma.pi_tilde(t,:) = ma.pi_T(t,:) - ma.pi_D(t,:)...                   % if the spv is bankrupt
                - sum(mi.r_li(t-1,:).*(mi.l_i(t-1,:)-mi.l_i_spvA(t-1,:))) ...   % equations are the same
                - sum(mi.r_li_spvA(t-1,:).*mi.l_i_spvA(t-1,:));...
 %              -
 %              sum(mi.r_li_spvA(t-1,div.bankrupt_spvA).*mi.l_i_spvA(t-1,div.bankrupt_spvA));
 %                  See above. Programm runs when code is commented out.
         end

% realized disposable income
            ma.y_d(t,:) = ma.omega(t,:)*ma.y(t,:) + ma.pi_D(t,:) + ...          % household disposable income consist of wage bill, distributed profits    
                          ma.r_spvA(t,:)*ma.spvA(t-1,:) + ...                         % interest income on spvA bonds                
                          ma.r_md(t,:)*ma.m_d_solvent(t-1,:);                          % interest income on deposits
            ma.eq_ch(t,:) = ma.pe_eb_solvent(t,:)- ma.pe_eb_solvent(t-1,:) +...         % capital gains on bank equities
                            ma.pe_ec_solvent(t,:)-ma.pe_ec_solvent(t-1,:);              % capital gains on coporate equities  
        
% household wealth
            ma.V_h(t,:) = ma.m_d(t-1,:) + ma.y_d(t,:) - ma.c(t,:) + ...
                        ma.pe_ec_solvent(t,:) + ma.pe_ec_insolvent(t,:) + ...
                        ma.pe_eb_solvent(t,:) + ma.pe_eb_insolvent(t,:) + ...
                        ma.spvA(t,:);              % change in household wealth equals the savings

% consumption share and profit rate            
ma.cons_share(t,:)          = ma.c(t,:)/ma.y(t,:);            
ma.prof_rate(t,:)           = ma.pi_T(t,:)/ma.K_solvent(t,:);
ma.y_d_share(t,:)           = ma.y_d(t,:)/ma.y(t,:);
ma.growth(t,:)              =(ma.y(t,:)-ma.y(t-1,:))/ma.y(t-1,:);
%% Compute further variables of interest
FurtherVars;                                                                    % determines fraction of bankrupt firms and banks                                                                      
            
 end

%% Checks (Consistency)
 Checks_balances;
%% Illustrations (figures)
% leverage structure of firm and banking sector
  stat.ER_k_mean = mean(mi.ER_k,2);                                           % returns the average equity ratio of banks
  stat.LR_mean   = mean(mi.LR_i,2);                                           % returns the average leverage ratio of firms
  stat.ER_k_5    = prctile(mi.ER_k,5,2);                                      % returns the five percent worst equity ratio
  stat.LR_5      = prctile(mi.LR_i,95,2);                                     % returns the five percent worst leverage ratio
%   Note: The last entry stands for the dimension, 
% "1" returns the mean/prctile of each column
% "2" returns the mean/prctile of each row