function [mi,div] = PSM_new(t,i,N_B,M,mi,~,div,LR_bar, par)
%% ------------------------------------------------------------------------
% Sort the banks arbitrarily in a vector          
            div.M_k = randperm(N_B);      
% Pick-up "M" banks randomly            
            div.M_k = div.M_k(1,1:M)';
           
%             'M' number of potential business partners
            div.r_lik = zeros(M,1);   % Initialize
            
%             Set a seed
%             rng(1,'twister')
            
            for ii = 1:M                                                     
                if mi.bankrupt_i(t,i)==0 
% Generate the interest rate offer
                        div.r_lik(ii,:) = div.r_lik(ii,:)...                    
                       + par.r_l0+ par.rho*(mi.LR_i(t,i)-LR_bar)...
                       -par.rho_ER*mi.ER_k(t-1,div.M_k(ii))+normrnd(par.mean_rl,par.sigma_rl); 
%                div.r_lik(ii,:) = div.r_lik(ii,:)...                    
%                        + par.r_l0+normrnd(par.mean_rl,par.sigma_rl); 
                else
% interest rate of insolvent firms
                 div.r_lik(ii,:) = 0;   
                end             
            end
% make sure that the stochasticity does not drive the interest rate below
% zero, i.e. determine lower threshold value
            div.r_lik(div.r_lik <0) = 0.001;                                
% Sort in ascending order            
           [div.r_sort, div.r_order] = sort(div.r_lik);                      
           div.M_ksort = div.M_k(div.r_order);              
% First element contains the bank offering the best credit condition (interest rate)           
           div.M_ksort = div.M_ksort';
           div.r_sort = div.r_sort';
%% Bank bankruptcies...  
            if mi.bankrupt_k(t,div.M_ksort(1,1)) == 1                
               counter = 2; 
% if a drawn bank goes bankrupt, take the second best bank,...
                while (mi.bankrupt_k(t,div.M_ksort(:,counter))==1) && (counter<M)
                    counter = counter+1;
                end
% If all banks in the sample are bankrupt, invest merely as much as it can finance by it's internal funds

                if (mi.bankrupt_k(t,div.M_ksort(:,counter))...
                        == mi.bankrupt_k(t,div.M_ksort(:,M)))...
                        && (mi.bankrupt_k(t,div.M_ksort(:,M))==1)
                   div.bestbank_i(t,i) = 0;
                   mi.transaction_cost(t,i) =0; 
                   mi.r_li(t,i) = mi.r_li(t-1,i); 
                else
                   div.bestbank_i(t,i) = div.M_ksort(1,counter);
                   mi.transaction_cost(t,i) = counter;
                   mi.r_li(t,i) =  div.r_sort(1,counter);
                end
            else 
            mi.r_li(t,i) =  div.r_sort(1,1);  
% Best Bank memory            
            div.bestbank_i(t,i) = div.M_ksort(1,1);                                  
            mi.transaction_cost(t,i) = 0;
           end
           
% For bankrupt firms
           if mi.bankrupt_i(t,i) == 1  
              div.bestbank_i(t,i) = 0;
           end
           
end

