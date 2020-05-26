%% Nur für grafische Illustrationen füe die Modellerklärung

% Scenario: 139
Time = 0:N-1;

fsize = 15;
range = 20;

figure('Visible','off')
plot(Time,MC.pi_spv(:,139),'Linewidth',2)
% hold on
% plot(res.g_Cup*100,'--k','Linewidth',1)
% hold on
% plot(res.g_Cbot*100,'--k','Linewidth',1)
% set(get(gca,'XLabel'),'Interpreter','latex','String','$$t$$','FontSize',fsize); 
% set(get(gca,'YLabel'),'Interpreter','latex','String','$$\pi^{spv}$$','FontSize',fsize,'Rotation',0); 
set(gca, 'Fontsize',fsize);
%  ylab = get(gca,'YLabel'); 
% set(ylab,'Position',get(ylab,'Position') - [0 .04 0])
xlim([0,range])
print -depsc SPVprofit;


% figure('Visible','off')
% plot(MC.ER_k(:,:,139),'Linewidth',2)
% set(gca, 'Fontsize',fsize);
% xlim([0,90])
% print -depsc ER;


figure('Visible','off')
hist(MC.ER_k(7,:,139))
set(gca, 'Fontsize',fsize);
xlim([0.11,0.22])
print -depsc ER_hist1;

figure('Visible','off')
hist(MC.ER_k(8,:,139))
set(gca, 'Fontsize',fsize);
xlim([0.11,0.22])
print -depsc ER_hist2;

figure('Visible','off')
plot(Time, MC.g(:,139),'Linewidth',2)
set(gca, 'Fontsize',fsize);
xlim([0,range])
print -depsc growth;

%% SPV Loans
div.vec_l = ones(N_C,1);
div.loans_spv = MC.l_i_spv(:,:,139)*div.vec_l;

figure('Visible','off')
plot(Time, div.loans_spv,'Linewidth',2)
set(gca, 'Fontsize',fsize);
xlim([0,range])
print -depsc loans_spv;
%% Transaction costs
div.vec_TC = ones(N_C,1);
div.transaction_costs = MC.transaction_cost(:,:,139)*div.vec_TC;

figure('Visible','off')
plot(Time, div.transaction_costs,'Linewidth',2)
set(gca, 'Fontsize',fsize);
xlim([0,range])
print -depsc transaction_costs;

%% Sonstige Grafiken
figure('Visible','off')
plot(Time, div.loans_spv,'Linewidth',2)
hold on
plot(Time, MC.bankrupt_firms(:,139),'Linewidth',2)
hold on
plot(Time, div.transaction_costs,'Linewidth',2)
hold on
plot(Time, MC.bankrupt_banks(:,139),'Linewidth',2)
set(gca, 'Fontsize',fsize);
xlim([0,range])
print -depsc loans_bankrupt;


figure('Visible','off')
plot(Time, MC.ER_k(:,:,139),'Color',[0.6,0.6,0.6],'Linewidth',2)
hold on
plot(Time, div.loans_spv/N_C,'Linewidth',2)
hold on
% plot(Time, MC.bankrupt_firms(:,139)/100,'Linewidth',2)
% set(gca, 'Fontsize',fsize);
xlim([0,range])
print -depsc ER_loans;
