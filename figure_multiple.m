%% Illustration

%% Define the framework

%results = struct('deposits',ma.m_d,'spv_profits',ma.pi_spvA,'bankrupt_banks',MC.bankrupt_banks,'credit',...
%    ma.l_c_solvent,'capital_stock',ma.K_solvent,'capital_stock_growth',ma.g,'consumption',ma.c,'gdp',ma.y,'gdp_growth',ma.growth);
% export series to excel
% 
% xlswrite('results_01_m_d.xlsx',ma.m_d);
% xlswrite('results_01_pi_spvA.xlsx',ma.pi_spvA);
% xlswrite('results_01_bankrupt_banks.xlsx',ma.bankrupt_banks);
% xlswrite('results_01_l_bank_solvent.xlsx',ma.l_c_solvent);
% xlswrite('results_01_K_solvent.xlsx',ma.K_solvent);
% xlswrite('results_01_i.xlsx',ma.i);
% xlswrite('results_01_c.xlsx',ma.c);
% xlswrite('results_01_y.xlsx',ma.y);
% xlswrite('results_01_growth.xlsx',ma.growth);
% 
% xlswrite('results_02_m_d.xlsx',ma.m_d);
% xlswrite('results_02_pi_spvA.xlsx',ma.pi_spvA);
% xlswrite('results_02_bankrupt_banks.xlsx',ma.bankrupt_banks);
% xlswrite('results_02_l_bank_solvent.xlsx',ma.l_c_solvent);
% xlswrite('results_02_K_solvent.xlsx',ma.K_solvent);
% xlswrite('results_02_i.xlsx',ma.i);
% xlswrite('results_02_c.xlsx',ma.c);
% xlswrite('results_02_y.xlsx',ma.y);
% xlswrite('results_02_growth.xlsx',ma.growth);
% 
% xlswrite('results_05_m_d.xlsx',ma.m_d);
% xlswrite('results_05_pi_spvA.xlsx',ma.pi_spvA);
% xlswrite('results_05_bankrupt_banks.xlsx',ma.bankrupt_banks);
% xlswrite('results_05_l_bank_solvent.xlsx',ma.l_c_solvent);
% xlswrite('results_05_K_solvent.xlsx',ma.K_solvent);
% xlswrite('results_05_i.xlsx',ma.i);
% xlswrite('results_05_c.xlsx',ma.c);
% xlswrite('results_05_y.xlsx',ma.y);
% xlswrite('results_05_growth.xlsx',ma.growth);
%  
% xlswrite('results_09_m_d.xlsx',ma.m_d);
% xlswrite('results_09_pi_spvA.xlsx',ma.pi_spvA);
% xlswrite('results_09_bankrupt_banks.xlsx',ma.bankrupt_banks);
% xlswrite('results_09_l_bank_solvent.xlsx',ma.l_c_solvent);
% xlswrite('results_09_K_solvent.xlsx',ma.K_solvent);
% xlswrite('results_09_i.xlsx',ma.i);
% xlswrite('results_09_c.xlsx',ma.c);
% xlswrite('results_09_y.xlsx',ma.y);
% xlswrite('results_09_growth.xlsx',ma.growth);


zz = zeros(30,3);
for i = 1:30
    zz(i)=i;
end
%results = struct('deposits',zz,'spv_profits',zz,'bankrupt_banks',zz,'credit',...
%    zz,'capital_stock',zz,'capital_stock_growth',zz,'consumption',zz,'gdp',zz,'gdp_growth',zz);
%import excel files
%results = table2struct(results1);

%import data as column vectors from results_gamma=085.xlsx
% results.deposits(:,1)       = results.m_d_01(:,1);
% results.spv_profits(:,1)    = results.pi_spvA_01(:,1);
% results.bankrupt_banks(:,1) = MC.bankrupt_banks_01(:,1);
% results.credit(:,1)         = results.l_c_solvent_01(:,1);
% results.capital_stock(:,1)  = results.K_solvent_01(:,1);
% results.capital_stock_growth(:,1) = results.g_01(:,1);
% results.consumption(:,1)    = results.c_01(:,1);
% results.gdp(:,1)            = results.y_01(:,1);
% results.gdp_growth(:,1)     = results.growth_01(:,1);
% 
% %import data from results_02.xlsx
% results.deposits(:,2)       = ma.m_d_02(:,1);
% results.spv_profits(:,2)    = ma.pi_spvA_02(:,1);
% results.bankrupt_banks(:,2) = MC.bankrupt_banks_02(:,1);
% results.credit(:,2)         = ma.l_c_solvent_02(:,1);
% results.capital_stock(:,2)  = ma.K_solvent_02(:,1);
% results.capital_stock_growth(:,2) = ma.g_02(:,1);
% results.consumption(:,2)    = ma.c_02(:,1);
% results.gdp(:,2)            = ma.y_02(:,1);
% results.gdp_growth(:,2)     = ma.growth_02(:,1);
% 
% %import data from results_05.xlsx
% results.deposits(:,3)       = ma.m_d_05(:,1);
% results.spv_profits(:,3)    = ma.pi_spvA_05(:,1);
% results.bankrupt_banks(:,3) = MC.bankrupt_banks_05(:,1);
% results.credit(:,3)         = ma.l_c_solvent_05(:,1);
% results.capital_stock(:,3)  = ma.K_solvent_05(:,1);
% results.capital_stock_growth(:,3) = ma.g_05(:,1);
% results.consumption(:,3)    = ma.c_05(:,1);
% results.gdp(:,3)            = ma.y_05(:,1);
% results.gdp_growth(:,3)     = ma.growth_05(:,1);

figure(1)
% plot results for securitization intensity = 0.1(blue), 0.2(green),
% 0.5(red)
subplot(3,3,1)
plot(zz,m_d_01,'b',zz,m_d_02,'g',zz,m_d_05,'r')
axis([15 30 50 250])
title('deposits (level)','FontSize',10)

subplot(3,3,2)
plot(zz,pi_spvA_01,'b',zz,pi_spvA_02,'g',zz,pi_spvA_05,'r')
axis([15 30 -0.1 0.1])
title('spv profits','FontSize',10)

subplot(3,3,3)
plot(zz,bankrupt_banks_01,'b',zz,bankrupt_banks_02,'g',zz,bankrupt_banks_05,'r')
axis([15 30 0 25])
title('bank insolvencies (%)','FontSize',10) 

subplot(3,3,4)
plot(zz,l_bank_solvent_01,'b',zz,l_bank_solvent_02,'g',zz,l_bank_solvent_05,'r')
axis([15 30 50 250])
title('credit (level)','FontSize',10)

subplot(3,3,5)
plot(zz,K_solvent_01,'b',zz,K_solvent_02,'g',zz,K_solvent_05,'r')
axis([15 30 300 450])
title('capital stock (level)','FontSize',10)

subplot(3,3,6)
plot(zz,i_01,'b',zz,i_02,'g',zz,i_05,'r')
axis([15 30 5 15])
title('net investment (level)','FontSize',10)

subplot(3,3,7)
plot(zz,c_01,'b',zz,c_02,'g',zz,c_05,'r')
axis([15 30 30 70])
title('consumption(level)','FontSize',10)

subplot(3,3,8)
plot(zz,y_01,'b',zz,y_02,'g',zz,y_05,'r')
axis([15 30 30 70])
title('GDP (level)','FontSize',10)

subplot(3,3,9)
plot(zz,growth_01,'b',zz,growth_02,'g',zz,growth_05,'r')
axis([15 30 -0.06 0.08])
title('GDP growth (%)','FontSize',10)

save results.mat
