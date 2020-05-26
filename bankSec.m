function [ mi,div] = bankSec( t,N_B, div, mi, ER_T, ~ )
%Optimize Banking Sector
   
    for k = 1: N_B
% For each bank "k", identify your business partners, i.e. the firms which 
% were drawn in the PSM
        div.node_k = find(div.bestbank_i(t,:)== k);                          
% Compute the bank-degrees        
        mi.bank_degrees(:,k) = numel(div.node_k); 
% Identify the respective credit conditions
        div.node_r_li = mi.r_li(t,div.node_k);                               
% ... and the respective size of the credit  
        div.mil_i_delta     = mi.l_i_delta(t,div.node_k);
% ...and the firms liabilities
        mi.r_times_l(t,k) = div.node_r_li*div.mil_i_delta';          
        

%% Compute bank's assets                                                       
        mi.l_s_ik(t,k) = mi.l_s_ik(t-1,k)+ sum(div.mil_i_delta);  
         if isempty(div.node_k)         % if a bank has no match in the PSM
            mi.l_s_ik(t,k) =  mi.l_s_ik(t-1,k); 
         end

    end
%% Bankruptcy

% mi.bankrupt_kk  = zeros(1,N_B); 
 
div.NCol    = find(mi.ER_k(t-1,:) < ER_T);
% if a bank was already bankrupt yesterday, it will stay bankrupt today
    mi.bankrupt_k(t,div.NCol) = 1;
%     mibankrupt_kk(:, mibankrupt_k(t-1,:)== 1 ) = 1;

% Banks equities
% This equation is really relevant for the development of the equity ratio
% Note that it is linked to the entire sector instead of the particular
% individual balance sheets and performances
 mi.p_ek(t,:)  =   mi.p_ek(t-1,:).*(1+ ((sum(mi.l_s_ik(t,:)) ...
     - sum(mi.l_s_ik(t-1,:)))./sum(mi.l_s_ik(t-1,:))) + normrnd(0,0.05));

% constant equity quantity
 mi.e_k(t,:)                 = mi.e_k(t-1,:);

% In case of bankruptcy, make sure that the individual bank equity has been
% frozen
mi.p_ek(t,div.NCol)   = (mi.l_s_ik(t,div.NCol) -  mi.m_d(t,div.NCol))...
    ./(mi.e_k(t,div.NCol));

% Here! Close balance sheet through md adjustment (residual)
% Deposits
 mi.m_d(t,:)   = mi.l_s_ik(t,:) -  mi.p_ek(t,:).* mi.e_k(t,:);
%% Equity Ratio
% The date "t" Equity Ratio
 mi.ER_k(t,:) = max(0,min(1,(mi.p_ek(t,:).*mi.e_k(t,:))./mi.l_s_ik(t,:))); 
% Implement it with a "t-1" index   
        mi.m_d(t,div.NCol)       =  mi.m_d(t-1,div.NCol);
        mi.p_ek(t,div.NCol)     =  mi.p_ek(t-1,div.NCol);
        mi.e_k(t,div.NCol)      =  mi.e_k(t-1,div.NCol);
        mi.ER_k(t,div.NCol)     =  mi.ER_k(t-1, div.NCol);
        



end

