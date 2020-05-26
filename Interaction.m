%% -------------------------------------------------------------------------
% Firm-Bank Interaction at the Micro Level
% -------------------------------------------------------------------------
% Note: The index "i" refers to firms, where "k" is the banking index.
%% Leverage Ratio
mi.LR_i(t,:)  =  mi.l_i(t-1,:)...                                           % Firm'S leverage ratio (as defined at the macro level)
                ./(mi.l_i(t-1,:)+ mi.p_ei(t-1,:).*mi.e_i(t-1,:));
[~, div.NCol] = find(mi.LR_i(t,:)  > LR_T);                                 % Determine insolvent firms: Find denotes a function and ~ is used for the empty case

%% Bankruptcy condition
% -------------------------------------------------------------------------
% Define a logical matrix (Boulsche string)
% if a firm is bankrupt, create a "1" for "bankrupt==true" in the binary 
% matrix and "0==not true" otherwise
% -------------------------------------------------------------------------
if ~isempty(div.NCol)                                                       % later div.NCol contains insolvent banks, here firms (non-empty)                                             
    mi.bankrupt_i(t,div.NCol) = 1;                                          % also consider firms as a flag
end
% The firm still remains bankrupt if it was already bankrupt yesterday
mi.bankrupt_i(t, mi.bankrupt_i(t-1,:)== 1 ) = 1;                            % the idea is that these firms do not participate in the market any longer

% Count the number of bankrupt firms in each period of time
% The appearing variable will be relevant for next period's macro-micro allocation.

div.non_operating_firms(t,:) = numel(find(mi.bankrupt_i(t,:) == 1));        % this is a time series containing the absolute number of insolvent firms
 
% Compute the firm shares of the market (concentration process)
 mi.Phi_i(t,:) =  1/(N_C-div.non_operating_firms(t,:));                           
if div.non_operating_firms(t,:) == N_C                                     
   mi.Phi_i(t,:) = 0;
end
% This loop runs for firms    
 for i = 1 : N_C   
                       
%% -------------------------------------------------------------
% Partner Selection Mechanism
% --------------------------------------------------------------                     
% calling the selection mechanism function
[mi,div] = PSM_new(t,i,N_B,M,mi,ma,div,LR_bar,par);
% The partner selection mechanism is a function of the time period t, the
% firm running index i, the search parameter M, the micro structure mi, the
% macro structure ma, the diverse structure div, the leverage ratio benchmark LR_bar,
% the parameter structure par. The output of the function is stored
% in the structures mi and div.
%% -------------------------------------------------------------
% Firm Sector
% --------------------------------------------------------------  
% calling the firmsector function
[ mi,ma,div ] = firmsector(t, i, mi, ma, div, d, par);
% The firmsector function is a function of the time period t, the firm
% running index i, the micro structure mi, the macro structure ma, the
% diverse structure div, the desired structure d and the parameter
% structure par. The output of the function is stored in the structures mi,
% ma and div.
% 
 end
%  Determine the amount of new loans the firm tries to take out from a 
%  commercial bank in order to finance the new amount of investment
mi.l_i_delta(t,:) = mi.l_i(t,:) - mi.l_i(t-1,:);  

%% -------------------------------------------------------------
% Banking Sector / SPV
% --------------------------------------------------------------
% calling the banksector function including securitization
if (t==t_spv) && (N_spv>0) 
[ mi,ma,div] = bankSec_SPV( t,t_spv,N_B, N_C, N_spv, N_spvA, div, mi, ma, ER_T, par);
% The bankSec_SPV function is a function of the time periods, the number of agents, the number of securitized loans, the div, mi, ma par structures 
% and the equity target ratio. There are essential elements in the diverse structure which store the securitized loans and debt costs depending on the firm flag.  
% In the moment of securitization SPV credits and lending rates are
% transferred from diverse to the micro structure.
    mi.l_i_spvA(t,div.firms_spvA)   = div.l_ki_spvA(:,div.firms_spvA);
    mi.r_li_spvA(t,div.firms_spvA)  = div.r_li_spvA(:,div.firms_spvA);
else                            
% banking sector for any other periods do not change SPV credits and
% lending rates and take output from bankSec function
    ma.spvA(t,:)        = ma.spvA(t-1,:);
    mi.l_i_spvA(t,:)    = mi.l_i_spvA(t-1,:);
    mi.r_li_spvA(t,:)   = mi.r_li_spvA(t-1,:);
    [ mi,div] = bankSec( t,N_B, div, mi, ER_T, par );
% In the case of non-securitisation, set the volume of the spvA tranche to the same as the previous period, 
% set the same for the securitised loans and their interest costs. All other variables stem from the procedure 
% of the banking sector without the addition of the special purpose vehicle.
end
% transfer bankrupt banks
if t+1<=T
       mi.bankrupt_k(t+1,:) =  mi.bankrupt_k(t,:);
end
%% Market Entry
% Sort Bank's equity ratios in a descending order in order to identify the
% the best and worst banks, where best and worst is measured in terms of
% size of the bank's equity ratio.
% Hence, the index of the "best" bank is stored as first element, ...
[ME.div.ER_sort(t,:), ME.div.ER_order(t,:)] = sort(mi.ER_k(t,:),'descend'); 

if t>par.zeta
% Identify a bank which is at least "zeta" periods bankrupt
    ME.div.mi.bankrupt_k = mi.bankrupt_k(t-par.zeta+1:t,:)';
% Bankrupt_k contains binary entries, whereby 1 = insolvent, 0 = solvent
    ME.count = ME.count';
% ME.count is a matrix, in which the individual columns correspond to the banks and the entries indicate how long a bank was insolvent.
    ME.count(:,t) = sum(ME.div.mi.bankrupt_k,2);
    ME.count = ME.count';
% Store the bank index in an auxiliary vector, if the bank is at least "zeta" periods bankrupt   
    ME.div.count_index = find(ME.count(t,:)>= par.zeta);                         
% If more than "eta" banks are contemporaneously bankrupt, choose just "eta"
% banks for recapitalisation. The others remain insolvent.
    if numel(ME.div.count_index)>=par.eta         
          ME.div.count_index =ME.div.count_index(:,randperm(numel(ME.div.count_index)));
          ME.div.count_index = ME.div.count_index(:,1:par.eta); 
          % Vector that contains the banks to be rescued/recapitalised
    end
      
    if ~isempty(ME.div.count_index)
% First: Switch the Boulsche number to "not-bankrupt" in order to allow the
% rescued bank to participate again at the Partner Selection Mechanism.
          mi.bankrupt_k(t,ME.div.count_index) = 0; 
            if t+1<= T
                mi.bankrupt_k(t+1,ME.div.count_index)   =   mi.bankrupt_k(t,ME.div.count_index);
            end
% Count the banks to be rescued (maximum number is eta)
          ME.div.numelcount = numel(ME.div.count_index);
% Identify the "best" banks which transfer a levy in terms of equities to 
% a special fund.
          ME.div.best_count = ME.div.ER_order(t,1:ME.div.numelcount);     
            
% The bank levies of the "best" banks are determined by a constant fraction 
% "kappa" of their equities
          ME.bestER(t,ME.div.best_count) = par.kappa*mi.e_k(t,ME.div.best_count);
% Carry out the equity transmission to the special fund and thereafter to
% the bankrupt bank
          mi.e_k(t,ME.div.best_count) = mi.e_k(t,ME.div.best_count) - ME.bestER(t,ME.div.best_count);
          mi.e_k(t,ME.div.count_index) = mi.e_k(t,ME.div.count_index) + ME.bestER(t,ME.div.best_count);
          mi.ER_k(t,ME.div.count_index) = (mi.p_ek(t,ME.div.count_index)...
              .*mi.e_k(t,ME.div.count_index) )./(mi.l_s_ik(t,ME.div.count_index)); 

     end

end


%% Lending rates
% Determine the macro interest rate on loans
ma.r_l(t,:)= (sum(mi.r_li(t,:) .* (mi.l_i(t,:)-mi.l_i_spvA(t,:)) ...
    ))/sum(mi.l_i(t,:)-mi.l_i_spvA(t,:));     

% Determine the interest rate on deposits (which is supposed to be less
% than the macro lending rate, i.e. theta_r is a multiplicative spread)
ma.r_md(t,:)= par.theta_r*(sum(mi.r_li(t,:) .* (mi.l_i(t,:)-mi.l_i_spvA(t,:))))...
    /sum(mi.l_i(t,:)-mi.l_i_spvA(t,:));     
% theta_r denotes the percentage of the lending rate that is paid for deposits 
% In case of securitization, identify the interest on SPV-bonds
    if sum(mi.l_i_spvA(t,:))>0
        ma.r_spvA(t,:) =  par.theta_spvA*(sum(mi.r_li_spvA(t,:).*mi.l_i_spvA(t,:))/sum(mi.l_i_spvA(t,:)));
    else
        ma.r_spvA(t,:) = 0;
    end
% theta_spvA is the percentage of the lending rate that is paid for the
% spvA bonds (tranche A): To make these bonds appealing for household
% theta_spvA > theta_r
        
%% Identify the bankrupt firms, in particular in the securitized loan portfolio
div.bankrupt_i = find(mi.bankrupt_i(t,:)==1);   
%   |no securitization| or |no bankruptcy in the SPV's loan portfolio|
if isempty(div.firms_spvA) || isempty(div.bankrupt_i)
        div.bankrupt_spvA = 0;
% div.bankrupt_spvA collects the firm flags of bankrupt firms in the
% securitized loan portfolio 
else
        div.bankrupt_spvA = zeros(1,numel(div.bankrupt_i));
    if numel(div.bankrupt_i) == 1 % if just one firm is insolvent,...
    % and if the respective firm is within the SPV's loan portfolio
    % -> Identify the respective firm
                div.bankrupt_spvA = find(div.bankrupt_i(:,1) == div.firms_spvA);
                % find those firms which are bankrupt and from which loans
                % have been securitized. it returns the index position in
                % div.firms_spvA, not the firm flag
                div.bankrupt_spvA = div.firms_spvA(:,div.bankrupt_spvA);
                % this stores the firm flag in div.bankrupt_spvA, as the
                % firm flag is the entry in div.firms_spvA(:,position)
    else % if more than one firm is bankrupt
                for s = 1 : numel(div.bankrupt_i) % for all loans of insolvent firms, if more than one firm is insolvent
% If it is in the SPV's loan portfolio, the rest is the same as above
                    if ~isempty(find(div.bankrupt_i(:,s) == div.firms_spvA, 1))
                        div.bankrupt_spvA(:,s) = find(div.bankrupt_i(:,s) == div.firms_spvA);
                        div.bankrupt_spvA(:,s) = div.firms_spvA(:,div.bankrupt_spvA(:,s));
                    else
                        div.bankrupt_spvA(:,s) = 0;  % "empty" else
                    end
                    div.bankrupt_spvA = div.bankrupt_spvA(:,div.bankrupt_spvA>0); 
                end
     end
end
% Finally, in order to make sure that the vector is not empty    
        if isempty(div.bankrupt_spvA)
            div.bankrupt_spvA = 0;
        end
%% SPV profits 
% are equal to debt service payments of all securitized loans minus those bankrupt (if there are bankrupt ones) minus interest payment to households (for tranche A) 
    if  div.bankrupt_spvA >0 
        ma.pi_spvA(t,:) = sum(mi.r_li_spvA(t,:).*mi.l_i_spvA(t,:))...               
        - sum(mi.r_li_spvA(t,div.bankrupt_spvA).*mi.l_i_spvA(t,div.bankrupt_spvA))...
        -ma.r_spvA(t,:)*ma.spvA(t,:);
    else
        ma.pi_spvA(t,:) = sum(mi.r_li_spvA(t,:).*mi.l_i_spvA(t,:))...     
        -ma.r_spvA(t,:)*ma.spvA(t,:); 
    end
%%  Liquiditation of the SPV 
% happens if the SPV-profits of the previous period are negative
    if ma.pi_spvA(t,:)<0
        ma.spvA(t,:)        = 0;
        mi.l_i_spvA(t,:)    = 0;
    end
    if  ma.pi_spvA(t-1,:)<0 
          ma.pi_spvA(t,:)   = 0;
          ma.spvA(t,:)      = 0;
        mi.l_i_spvA(t,:)    = 0;
    end
       
    if ma.pi_spvA(t,:) <0
%  Exclude loans from insolvent firms
div.firms_spvA_solvent=div.firms_spvA;
div.firms_spvA_solvent(ismember(div.firms_spvA,div.bankrupt_spvA))=[];
        
% Identify the banks of the SPV loan protfolio        
        div.bank_spv                = div.bestbank_i(t_spv,div.firms_spvA_solvent);
% Avoid double counting
        div.unique_banks            = unique(div.bank_spv);
% Theoretically, the same bank could be bestbank for more than one firm
        div.count_banks             = hist(div.bank_spv,div.unique_banks);
        div.indexToRepeatedValue    = (div.count_banks~=1);
        div.bank_spv                = unique(div.bank_spv);


%% ------------------------ SPV Bankruptcy ------------------------------------

% Here, differentiate between the losses of the banks and the investors.
% Distribute the losses via the share gamma to both, households and banks.

% Rückbuchung der Kredite
mi.l_s_ik(t,div.bank_spv) =  mi.l_s_ik(t,div.bank_spv)+ div.l_ki_spvA_NB(t_spv,div.bank_spv)./5;

% gamma regelt ob den zurückgebuchten Krediten Depositen oder Bankaktien gegenüberstehen 
mi.p_ek(t,div.bank_spv)   = par.gamma*mi.p_ek(t,div.bank_spv).*(1+ ((sum(mi.l_s_ik(t,:)) - sum(mi.l_s_ik(t-1,:)))./sum(mi.l_s_ik(t-1,:))))...
    + (1-par.gamma)*(mi.p_ek(t,div.bank_spv).*mi.e_k(t,div.bank_spv) + div.l_ki_spvA_NB(t_spv,div.bank_spv))./mi.e_k(t,div.bank_spv);
%mi.p_ek(t,:)  =   mi.p_ek(t-1,:).*(1+ ((sum(mi.l_s_ik(t,:)) - sum(mi.l_s_ik(t-1,:)))./sum(mi.l_s_ik(t-1,:))));


% Adjust Deposits 
mi.m_d(t,:)   = mi.l_s_ik(t,:) -  mi.p_ek(t,:).* mi.e_k(t,:);
 
% Adjust equity ratios
mi.ER_k(t,div.bank_spv) = max(0,min(1,(mi.p_ek(t,div.bank_spv).*mi.e_k(t,div.bank_spv)./(mi.l_s_ik(t,div.bank_spv)))));

% Make sure that the banks which were already bankrupt are not affected
%         mi.m_d(t,div.NCol)       =  mi.m_d(t-1,div.NCol);
%         mi.p_ek(t,div.NCol)     =  mi.p_ek(t-1,div.NCol);
%         mi.e_k(t,div.NCol)      =  mi.e_k(t-1,div.NCol);
%         mi.ER_k(t,div.NCol)     =  mi.ER_k(t-1, div.NCol);
% Check balance Sheet
% check.bank_arbitrary = [mi.l_s_ik(t,div.bank_spv(1,1)),...
%     mi.p_ek(t,div.bank_spv(1,1))*mi.e_k(t,div.bank_spv(1,1)) + mi.m_d(t,div.bank_spv(1,1))];
    end

% Set spv variables equal to 0 after liquidation        
    if (ma.spvA(t-1,:) == 0) && (t>t_spv)
        ma.spvA(t,:) = ma.spvA(t-1,:);
        mi.l_i_spvA(t,:) = mi.l_i_spvA(t-1,:);
    end
    
%% Determine level of loans after interaction and bankruptcies
% Firm sector
    ma.l_c_insolvent(t,:) = sum(mi.l_i(t, mi.bankrupt_i(t,:) ==1));
% Banking sector
    ma.l_bank_insolvent(t,:) = sum(mi.l_s_ik(t,mi.bankrupt_k(t,:) == 1));
    ma.l_bank_solvent(t,:)   = sum(mi.l_s_ik(t,mi.bankrupt_k(t,:) ==0));

% Credit demand after interaction
    ma.l_c_solvent(t,:) = sum(mi.l_i(t,:)); % Note that these are all corporate loans, not only those of solvent firms
