function [ mi,ma,div ] = firmsector(t, i, mi, ma, div, d, par)
%New Attempt to optimize the firm sector
% For all solvent fimrs, compute...
           if  mi.bankrupt_i(t,i) == 0 
               mi.i_i(t,i) = mi.Phi_i(t,:)*d.i(t,:); 
               mi.K_i(t,i) = mi.i_i(t,i) + (1-par.delta)*mi.K_i(t-1,i);
               
               mi.pi_i(t,i)  = mi.Phi_i(t,:)*d.pi_T(t,:); 
               mi.pi_iD(t,i) = par.xi_c*mi.pi_i(t,i);                                       
               mi.pi_tilde_i(t,i) = mi.pi_i(t,i) - mi.pi_iD(t,i) ...
                   -mi.r_li(t-1,i)*(mi.l_i(t-1,i)-mi.l_i_spvA(t-1,i))...
                   - mi.r_li_spvA(t-1,i)*mi.l_i_spvA(t-1,i);
%% Searchcosts/Transactioncosts    
                if mi.transaction_cost(t,i) == 0
                    mi.l_i(t,i)  = mi.l_i(t-1,i) + mi.i_i(t,i) -mi.pi_tilde_i(t,i);    
                else % Searchcosts
                    mi.l_i(t,i)  = mi.l_i(t-1,i) + mi.i_i(t,i) -mi.pi_tilde_i(t,i);
                    mi.l_i(t,i)  = mi.l_i(t,i)...
                        - (mi.transaction_cost(t,i)/100)*par.searchcost*mi.l_i(t,i); 
%                   Adjust firms investment by the searchcosts
                    mi.i_i(t,i) = mi.l_i(t,i) - mi.l_i(t-1,i)+mi.pi_tilde_i(t,i);
                    mi.K_i(t,i)= mi.i_i(t,i) + (1-par.delta)*mi.K_i(t-1,i);
                end
                if t <=2
                     mi.var_pei(t,i) = par.sig0 + par.sig1*((mi.p_ei(t-1,i) - mi.p_ei(1,1)))^2 ...
                        + par.sig2 * mi.var_pei(t-1,i) ;   
                    mi.p_ei(t,i) =  abs(normrnd(mi.p_ei(1,1),par.sig0));
                else
                     mi.var_pei(t,i) = par.sig0 + par.sig1*(log(mi.p_ei(t-1,i) / mi.p_ei(t-2,i)))^2 ...
                        + par.sig2 * mi.var_pei(t-1,i) ; 
                     mi.sigma_pei(t,i) = sqrt(mi.var_pei(t,i));

                    mi.p_ei(t,i) = exp( ...
                    log(mi.p_ei(t-1,i))+mi.sigma_pei(t,i)*normrnd(0,1));
                end
                
    
    mi.e_i(t,i) = (mi.K_i(t,i) - mi.l_i(t,i))/mi.p_ei(t,i); 
%% Transaction costs after internal financing
 mi.i_i(t,i) = mi.l_i(t,i)-mi.l_i(t-1,i) + mi.pi_tilde_i(t,i);
 mi.K_i(t,i) = mi.i_i(t,i) + (1-par.delta)*mi.K_i(t-1,i);
 %% Neu!!! Schließt die Corporate Balance
 mi.e_i(t,i) = (mi.K_i(t,i) - mi.l_i(t,i))/mi.p_ei(t,i); 
           else  % insolvent firms (no investment)
               mi.i_i(t,i) = 0;
               mi.K_i(t,i) = mi.K_i(t-1,i);
               mi.pi_i(t,i)  = 0; 
               mi.pi_iD(t,i) = 0;                                              
               mi.pi_tilde_i(t,i) = 0;%mi.pi_i(t,i) - mi.pi_iD(t,i) -mi.r_li(t-1,i)*mi.l_i(t-1,i);
               mi.l_i(t,i)  = mi.l_i(t-1,i);
               mi.e_i(t,i) = (mi.K_i(t,i) - mi.l_i(t,i))/mi.p_ei(t,i); 
               mi.p_ei(t,i) = mi.p_ei(t-1,i);
           end
          
           %%
% Adjust investment of the firms who found no bank (partner) in the PSM due
% to bank bankruptcies
            mi.l_rest(t,i) = 0;
            if (div.bestbank_i(t,i) ==0) && (mi.bankrupt_i(t,i) == 0)
% In such case, investment are financed by retained earnings                
                mi.i_i(t,i) = mi.pi_tilde_i(t,i);
                mi.K_i(t,i) = mi.i_i(t,i) + (1-par.delta)*mi.K_i(t-1,i);
                mi.l_i(t,i)  = mi.l_i(t-1,i); 
                mi.l_rest(t,i)      = mi.l_i(t,i);
                %%
                if t <=2
                     mi.var_pei(t,i) = par.sig0 + par.sig1*((mi.p_ei(t-1,i) - mi.p_ei(1,1)))^2 ...
                        + par.sig2 * mi.var_pei(t-1,i) ;   
                    mi.p_ei(t,i) =  abs(normrnd(1,0.01));
                else
                     mi.var_pei(t,i) = par.sig0 + par.sig1*(log(mi.p_ei(t-1,i) / mi.p_ei(t-2,i)))^2 ...
                        + par.sig2 * mi.var_pei(t-1,i) ; 
                     mi.sigma_pei(t,i) = sqrt(mi.var_pei(t,i));

                    mi.p_ei(t,i) = exp( ...
                    log(mi.p_ei(t-1,i))+mi.sigma_pei(t,i)*normrnd(0,1));
                end
                %% Neu!!! Schließt die Corporate Balance
 mi.i_i(t,i) = mi.l_i(t,i)-mi.l_i(t-1,i) + mi.pi_tilde_i(t,i);
 mi.K_i(t,i) = mi.i_i(t,i) + (1-par.delta)*mi.K_i(t-1,i);
 mi.e_i(t,i) = (mi.K_i(t,i) - mi.l_i(t,i))/mi.p_ei(t,i); 
            elseif (div.bestbank_i(t,i) ==0) && (mi.bankrupt_i(t,i) == 1)
                mi.i_i(t,i) = 0;
                mi.K_i(t,i) = mi.K_i(t-1,i);
                mi.pi_i(t,i)  = 0; 
                mi.pi_iD(t,i) = 0;                                              
                mi.pi_tilde_i(t,i) = 0; %mi.pi_i(t,i) - mi.pi_iD(t,i) -mi.r_li(t-1,i);
                mi.l_i(t,i)  = mi.l_i(t-1,i); 
                mi.p_ei(t,i) = mi.p_ei(t-1,i); 
                mi.e_i(t,i) = (mi.K_i(t,i) - mi.l_i(t,i))/mi.p_ei(t,i); 
            end

end

