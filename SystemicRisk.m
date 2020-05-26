%% Keep the house clean
clear;
clc;
close all;


%% Cases
Noofcases =3;

VaR.percy = 5; % one side


for j = 1: Noofcases
%% load data (from simulations)
if      j == 1
    load Workspaces/WSs/21112017/MC1_workspace1000.mat 
    
    %% Compute mean of all equity ratio over all banks (mean of each column)
VaR.mean_time_ER_k_1 = zeros(Ncol,N_B);
% fill the matrix with 
% -replications in rows
% -banks in columns
for i = 1:Ncol
    VaR.mean_time_ER_k_1(i,:) = nanmean(MC.ER_k(:,:,i),1);
end

VaR.mean_time_bank_ER_k_1 = nanmean(VaR.mean_time_ER_k_1,2);

%% Compute the x% percentile
% sort in a ascending order (defualt); the first entries are the lowest 
% mean equity ratio and thus the worst.
VaR.percy_entry = (VaR.percy/100)*Ncol;

   [dif.VaRsort,dif.VaRorder] = sort(VaR.mean_time_bank_ER_k_1);
   
VaR.VaR_5percy_1      = dif.VaRsort(VaR.percy_entry,:); 
VaR.VaR_5difference_1 = dif.VaRsort(VaR.percy_entry,:) -dif.VaRsort(1,:);
VaR.VaR_5mean_1       = nanmean(dif.VaRsort(1:VaR.percy_entry,:));


elseif  j == 2
    load Workspaces/WSs/21112017/MC2_workspace1000.mat 
    
    %% Compute mean of all equity ratio over all banks (mean of each column)
VaR.mean_time_ER_k_5 = zeros(Ncol,N_B);
% fill the matrix with 
% -replications in rows
% -banks in columns
for i = 1:Ncol
    VaR.mean_time_ER_k_5(i,:) = nanmean(MC.ER_k(:,:,i),1);
end

VaR.mean_time_bank_ER_k_5 = nanmean(VaR.mean_time_ER_k_5,2);

%% Compute the x% percentile

   [dif.VaRsort,dif.VaRorder] = sort(VaR.mean_time_bank_ER_k_5);
   
VaR.VaR_5percy_5      = dif.VaRsort(VaR.percy_entry,:); 
VaR.VaR_5difference_5 = dif.VaRsort(VaR.percy_entry,:) -dif.VaRsort(1,:);
VaR.VaR_5mean_5       = nanmean(dif.VaRsort(1:VaR.percy_entry,:));
else    
    load Workspaces/WSs/21112017/MC3_workspace1000.mat 
    %% Compute mean of all equity ratio over all banks (mean of each column)
VaR.mean_time_ER_k_10 = zeros(Ncol,N_B);
% fill the matrix with 
% -replications in rows
% -banks in columns
for i = 1:Ncol
    VaR.mean_time_ER_k_10(i,:) = nanmean(MC.ER_k(:,:,i),1);
end

VaR.mean_time_bank_ER_k_10 = nanmean(VaR.mean_time_ER_k_10,2);

%% Compute the x% percentile

   [dif.VaRsort,dif.VaRorder] = sort(VaR.mean_time_bank_ER_k_10);
   
VaR.VaR_5percy_10      = dif.VaRsort(VaR.percy_entry,:); 
VaR.VaR_5difference_10 = dif.VaRsort(VaR.percy_entry,:) -dif.VaRsort(1,:);
VaR.VaR_5mean_10       = nanmean(dif.VaRsort(1:VaR.percy_entry,:));
    
    
     save VaR_measures VaR;

end
   

 
end