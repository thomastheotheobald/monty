function [ mi,ma,div] = SPV_new( t, t_spv,N_B, N_C, N_spv,N_spvA, mi,ma, div)
%Optimized Version
%% --------------------------------------------------------
% Securitization
% ---------------------------------------------------------

%% Randomly Sample loans
% ... and make sure that the sample excludes lemons
    div.M_spv       = 1:N_C;
    div.M_not_bankrupt = div.M_spv(:,div.bestbank_i(t,:)>0);

if (isempty(div.M_not_bankrupt)) || (numel(div.M_not_bankrupt)< N_spv)  ...
        || (N_spv==0)    % if all in the sample are bankrupt
%     ma.spvA = 0;
    div.r_li_spvA = 0;
    div.l_ki_spvA = 0;
    div.firms_spvA = 0;
else
% determine randomly the loans which are subject of the securitization    
    div.M_not_bankrupt = div.M_not_bankrupt(:,randperm(numel(div.M_not_bankrupt)));
    div.M_spv          = div.M_not_bankrupt(:,1:N_spv);

%  determine the loans which are subject of the securitization
    div.l_c_spv     = mi.l_i(t,div.M_spv);  
% In order to assess the qualitiy of the loans, pick-up the respective firm
% leverage ratios    
    div.LR_ispv     = mi.LR_i(t,div.M_spv);                                 
%%  Sort from best to worse (ascending order)
    if t == 1
        div.LR_order = div.M_spv;                                           % random distribution in period 1
        div.LR_sort  = mi.LR_i(t,div.LR_order);
    else
        [div.LR_sort, div.LR_order] = sort(div.LR_ispv);                    
    end 
%  Sort the respective loans according to the leverage ratios  
        div.M_order = div.M_spv(div.LR_order);
        mi.l_hat = mi.l_i(t,div.M_order);                               
%%    
%  Determine the A-tranche
    if  N_spv          == 0
        ma.spvA(t,:)         = 0;
    else
        ma.spvA(t,:)      =  sum(mi.l_hat(:,1:N_spvA));                           % SPV-A Tranche (macro)

        div.firms_spvA = div.M_order(1:N_spvA);                             % Identify the firms
%         div.l_ki_spvA(:,div.firms_spvA)  = mi.l_i_delta(t,div.firms_spvA);  % ... and loans
        div.l_ki_spvA(:,div.firms_spvA)  = mi.l_i(t,div.firms_spvA);  % ... and loans

        div.r_li_spvA(:,div.firms_spvA)  = mi.r_li(t,div.firms_spvA);       % ... and interest rates
% ... and determine the respective commercial bank behind that loans        
        div.bank_spvA  = div.bestbank_i(t, div.firms_spvA);                  
% Make sure that the banks just count once and allocate the credits
allocate_loans = zeros(1,N_B);
    for bb = div.bank_spvA
        BankFirmNetwork_Logical = div.bestbank_i(t_spv,:)==bb;
%         disp(BankFirmNetwork_Logical)
%         BankFirm_Loans          = find(BankFirmNetwork_Logical==1); 
%         allocate_loans(:,bb) = allocate_loans(:,bb)+...
%             sum(div.l_ki_spvA(:,div.bank_spvA==bb));  
        allocate_loans(:,bb) = allocate_loans(:,bb)+...
            mi.l_i(t_spv,:)*BankFirmNetwork_Logical';
    end
    
    
%% Create Vector with banks
        div.unique_banks = unique(div.bank_spvA);
        div.count_banks = hist(div.bank_spvA,div.unique_banks);
        div.indexToRepeatedValue = (div.count_banks~=1);
        
        div.bank_spvA = div.unique_banks;
% Take out the securitized A-Tranche
        div.l_ki_spvA_NB(t,:) = allocate_loans;    
        
        mi.l_s_ik(t,div.bank_spvA) = max(mi.l_s_ik(t,div.bank_spvA) - div.l_ki_spvA_NB(t,div.bank_spvA)./5,0);                             

   end
end    

end

