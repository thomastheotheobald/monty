
%bank balance (solvent banks)
check.BB_solvent = [ma.l_bank_solvent , ma.pe_eb_solvent+ ma.m_d_solvent];  
%bank balance (insolvent banks)
check.BB_insolvent = [ma.l_bank_insolvent, ma.pe_eb_insolvent + ma.m_d_insolvent];
%corporate balance (solvent corporates)
check.CB_solvent = [sum(mi.K_i,2), ma.pe_ec_solvent + ma.l_c_solvent];
%corporate balance (insolvent corporates)
check.CB_insolvent = [ma.K_insolvent, ma.pe_ec_insolvent + ma.l_c_insolvent]; 
%savings and investment
check.SI = [ma.y_d  - ma.c + ma.pi_tilde, ma.i];
%wealth without capital gains
check.V_check = [ma.V_h, ma.m_d_solvent + ma.spvA + ma.p_ec.*ma.e_c + ma.p_eb.*ma.e_b];
%disposable income
check.disp_check = [ma.y_d - (ma.pe_eb_solvent-ma.pe_eb_solvent) - (ma.pe_ec_solvent-ma.pe_ec_solvent)+ma.pi_tilde,ma.y];
e_c_growth = zeros(T,1);
K_growth   = zeros(T,1);
Loan_growth       = zeros(T,1); 
% L_sik = zeros(T,1);

% for z =2 :T
%     e_c_growth(z) = (ma.p_ec(z) - ma.p_ec(z-1))/ma.p_ec(z-1);
%     K_growth(z) = (sum(mi.K_i(z,mi.bankrupt_i(z,:)==0)) - ...
%         sum(mi.K_i(z-1,mi.bankrupt_i(z-1,:)==0)))/sum(mi.K_i(z-1,mi.bankrupt_i(z-1,:)==0));
% %     Loan_growth(z) = (ma.l_c_solvent(z) - ma.l_c_solvent(z-1))/ma.l_c_solvent(z-1);
%     Loan_growth(z)    = (sum(mi.l_s_ik(z,:)) - sum(mi.l_s_ik(z-1,:)))/sum(mi.l_s_ik(z-1,:));
% end
L_sik = sum(mi.l_s_ik,2);
L_i   = sum(mi.l_i,2); 
L_i_delta = sum(mi.l_i_delta,2);
% check.GC = [K_growth , e_c_growth, Loan_growth];

% L_check = [ma.l_bank_solvent, ma.l_c_solvent, BB_solvent(:,2)];
% L_check = [ma.l_bank_solvent, ma.l_c_solvent, L_sik];
check.L_check = [L_sik, L_i,ma.l_c_solvent];


check.K_check = [ma.K, ma.K_solvent, ma.K_insolvent];
 

check.Pi_T_check = [ma.pi_T, ma.pi_T_solvent];

check.R_check = [ma.r_l, ma.r_spvA ,ma.r_md ];
% figure
% plot(mi.K_i(70:end,:))

mi.CB = zeros(T,2);
% for o = 1:N_C
%    mi.CB(:,:,o) = [mi.K_i(:,o), mi.l_i(:,o) + mi.p_ei(:,o)*mi.e_i]; 
% end
% miCB
% AAAA = miCB(:,:,1);
%% Individuelle Bank-Checks

% partiell
check.bank_k_asset = mi.l_s_ik;
check.bank_k_liability = mi.p_ek.*mi.e_k + mi.m_d;

% % check.bank_k = zeros(T,N_B);
% % for o = 1:N_B
% %     check.bank_k(:,o) = check.bank_k_equity(:,o) - check.bank_k_liability(:,o);
% % end
check.bank_k = check.bank_k_asset - check.bank_k_liability;

check.firm_i = mi.K_i - mi.l_i- mi.p_ei.*mi.e_i;




%% 
div.g_l = zeros(numel(ma.l_c_solvent(:,1)),1);
div.g_K = zeros(numel(ma.l_c_solvent(:,1)),1);


for o = 2:length(ma.l_c_solvent(:,1))
div.g_l(o) = div.g_l(o)+(ma.l_c_solvent(o,:) - ma.l_c_solvent(o-1,:))/ma.l_c_solvent(o-1,:);
div.g_K(o) = div.g_K(o)+(ma.K(o,:) - ma.K(o-1,:))/ma.K(o-1,:);

end
check. g_KL = [div.g_l ,div.g_K];

check.CB_alternative = [ma.K, ma.l_c_solvent+ma.pi_tilde];


% Interest rates

check.r = [ma.r_l, ma.r_md, ma.r_spvA];