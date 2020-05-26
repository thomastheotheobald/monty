%% Für illustrative Zwecke
% Identifiziere die SPV Pleite Periode
if (ma.pi_spvA(t-1,:) ~= 0) && (ma.pi_spvA(t,:) == 0)
   MC.spv_bankrupt = t;  
end
    

%% Diese hier vielleicht ausgliedern  
            ma.bankrupt_firms(t,:) = ma.bankrupt_firms(t,:)...              % Determine fraction of bankrupt firms
                +(numel(find(mi.bankrupt_i(t,:) ==1))/N_C)*100;
            ma.bankrupt_banks(t,:) = ma.bankrupt_banks(t,:)...              % Determine fraction of bankrupt banks
                +(numel(find(mi.bankrupt_k(t,:) ==1))/N_B)*100;