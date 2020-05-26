%% Keep the house clean


%% Load the Benchmark Scenario and create the figures
% load Workspaces/BenchmarkWS/16022017_v3.mat;


%% Define time index (x-axis)
  TimeIndex = 1:N;
  TimeIndex = (TimeIndex-1)';
  
%   horizone
  per = 29;
%   per = per;
  fontsz = 13;
 %% figure I-S 
%  put the check in the corresponding folder
check.SI = MC.SI;
 
    figure('Visible','off')
h1 = plot(TimeIndex(1:per+1),check.SI(1:per+1,1),'color',[0 0.447 0.7410], 'Linewidth', 2);
hold on
h2 = plot(TimeIndex(1:per+1),check.SI(1:per+1,2), 'Linewidth', 1);
%   ylim([0.015,0.05]);
% set(get(gca,'YLabel'),'Interpreter','latex','String','$$g_{t}$$'...
%         ,'Units', 'Normalized', 'Position', [-0.1, 0.46, 0],'FontSize',18); 
set(get(gca,'YLabel'),'Interpreter','latex','String','Investment, Saving','FontSize',fontsz);
set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',fontsz);
% set(get(gca,'YLabel'),'Rotation',0)
  set(gca, 'Fontsize',13)
hl=legend([h1 h2],{'Saving','Investment'});
 set(hl, 'interpreter', 'latex', 'Location', 'northwest')
 legend('boxoff')
  print -depsc figures/Checks/SI;
  
  %% Firm-Sector (solvent firms)
  clearvars check;
  
  check.solvent_firms= find(MC.bankrupt_i(end,:)==0);
  check.CB_assets_solvent = mi.K_i(cut:end,check.solvent_firms);
  check.CB_liabilities = mi.l_i(cut:end,check.solvent_firms)...
      + (mi.p_ei(cut:end,check.solvent_firms).*mi.e_i);

      figure('Visible','off')
h1 = plot(TimeIndex(1:per+1),sum(check.CB_assets_solvent(1:per+1,1),2),'color',[0 0.447 0.7410], 'Linewidth', 2);
hold on
h2 = plot(TimeIndex(1:per+1),sum(check.CB_liabilities(1:per+1,1),2), 'Linewidth', 1);
%   ylim([0.015,0.05]);
% set(get(gca,'YLabel'),'Interpreter','latex','String','$$g_{t}$$'...
%         ,'Units', 'Normalized', 'Position', [-0.1, 0.46, 0],'FontSize',18); 
set(get(gca,'YLabel'),'Interpreter','latex','String','Balance Sheets of Solvent Firms','FontSize',fontsz);
set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',fontsz);
% set(get(gca,'YLabel'),'Rotation',0)
  set(gca, 'Fontsize',13)
hl=legend([h1 h2],{'Assets$$^{i}$$','Liabilities$$^{i}$$'});
 set(hl, 'interpreter', 'latex', 'Location', 'northwest')
 legend('boxoff')
  print -depsc figures/Checks/CB_solvent;
  
%% Firm-Sector (insolvent firms)

  check.insolvent_firms             = find(MC.bankrupt_i(end,:)==1);
  check.CB_assets_insolvent          = mi.K_i(cut:end,check.insolvent_firms);
  check.CB_liabilities_insolvent    = mi.l_i(cut:end,check.insolvent_firms)...
      + (mi.p_ei(cut:end,check.insolvent_firms)*mi.e_i);  
  
  figure('Visible','off')
h1 = plot(TimeIndex(1:per+1),sum(check.CB_assets_insolvent(1:per+1,1),2),'color',[0 0.447 0.7410], 'Linewidth', 2);
hold on
h2 = plot(TimeIndex(1:per+1),sum(check.CB_liabilities_insolvent(1:per+1,1),2), 'Linewidth', 1);
set(get(gca,'YLabel'),'Interpreter','latex','String','Balance Sheets of Insolvent Firms','FontSize',fontsz);
set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',fontsz);
  set(gca, 'Fontsize',13)
hl=legend([h1 h2],{'Assets$$^{i}$$','Liabilities$$^{i}$$'});
 set(hl, 'interpreter', 'latex', 'Location', 'southeast')
 legend('boxoff')
  print -depsc figures/Checks/CB_insolvent;
  
  %% Banking sector ((in- &)solvent Banks)
  clearvars div;

%   dimension
  [xb,yb]=size(MC.bankrupt_k);
%   initialize
  check.solvent_banks               = zeros(xb,yb);
  check.BB_assets_solvent           = zeros(xb,yb);
  check.BB_deposits_solvent         = zeros(xb,yb);
  check.BB_equities_solvent         = zeros(xb,yb);
  check.BB_liabilities_solvent      = zeros(xb,yb);
 
  check.BB_assets_insolvent         = zeros(xb,yb);
  check.BB_deposits_insolvent       = zeros(xb,yb);
  check.BB_equities_insolvent       = zeros(xb,yb);
  check.BB_liabilities_insolvent    = zeros(xb,yb);


%   Wash out transient periods (length of 'cut')
  MC.l_s_ik = mi.l_s_ik(cut:end,:);
  MC.m_d    = mi.m_d(cut:end,:);
  MC.e_k    = mi.e_k(cut:end,:);
  MC.p_ek   = mi.p_ek(cut:end,:);

% Compute the cross-product
  MC.pe_ek = MC.p_ek.*MC.e_k;
  
  for t = 1:xb
%% loans solvent banks

[div.r, div.c] = find(MC.bankrupt_k(t,:)==1);

% Identify the banks of the SPV loan protfolio        
div.banks_insolvent = unique(div.c);
% Identify the solvent banks
div.banks   = 1:N_B;
div.banks_solvent   = div.banks(~ismember(div.banks,div.banks(div.banks_insolvent)));
div.banks_solvent_boul  = zeros(1,yb);
    div.banks_solvent_boul(:,div.banks_solvent) = 1;
check.BB_assets_solvent(t,:) = check.BB_assets_solvent(t,:)...
    + MC.l_s_ik(t,:).*div.banks_solvent_boul;

% Deposits
  check.BB_deposits_solvent(t,:) =  check.BB_deposits_solvent(t,:)...
  + MC.m_d(t,:).*div.banks_solvent_boul;
% Equities (cross-product)
  check.BB_equities_solvent(t,:) =  check.BB_equities_solvent(t,:)...
  + MC.pe_ek(t,:).*div.banks_solvent_boul;


%% loans insolvent banks
div.banks_insolvent_boul  = zeros(1,yb);
    div.banks_insolvent_boul(:,div.banks_insolvent) = 1;

check.BB_assets_insolvent(t,:) = check.BB_assets_insolvent(t,:)...
    + MC.l_s_ik(t,:).*div.banks_insolvent_boul;

% Deposits
  check.BB_deposits_insolvent(t,:) =  check.BB_deposits_insolvent(t,:)...
  + MC.m_d(t,:).*div.banks_insolvent_boul;
% Equities (cross-product)
  check.BB_equities_insolvent(t,:) =  check.BB_equities_insolvent(t,:)...
  + MC.pe_ek(t,:).*div.banks_insolvent_boul;
  end

check.BB_liabilities_solvent = check.BB_equities_solvent...
    +check.BB_deposits_solvent ;
check.BB_liabilities_solvent_aggregate = sum(check.BB_liabilities_solvent,2);
  
check.BB_liabilities_insolvent = check.BB_equities_insolvent...
    +check.BB_deposits_insolvent;
check.BB_liabilities_insolvent_aggregate = sum(check.BB_liabilities_insolvent,2); 

check.BB_assets_solvent_aggregate = sum(check.BB_assets_solvent,2);
check.BB_assets_insolvent_aggregate = sum(check.BB_assets_insolvent,2);
check.BB_deposits_solvent_aggregate = sum(check.BB_deposits_solvent,2);



%% Figure

      figure('Visible','off')
h1 = plot(TimeIndex(1:per+1),check.BB_assets_solvent_aggregate(1:per+1,1),'color',[0 0.447 0.7410], 'Linewidth', 2);
hold on
h2 = plot(TimeIndex(1:per+1),check.BB_liabilities_solvent_aggregate(1:per+1,1), 'Linewidth', 1);
 
set(get(gca,'YLabel'),'Interpreter','latex','String','Balance Sheets of Solvent Banks','FontSize',fontsz);
set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',fontsz);
% set(get(gca,'YLabel'),'Rotation',0)
  set(gca, 'Fontsize',13)
hl=legend([h1 h2],{'Assets$$^{k}$$','Liabilities$$^{k}$$'});
 set(hl, 'interpreter', 'latex', 'Location', 'northwest')
 legend('boxoff')
  print -depsc figures/Checks/BB_solvent;
  
  
  figure('Visible','off')
h1 = plot(TimeIndex(1:per+1),check.BB_assets_insolvent_aggregate(1:per+1,1),'color',[0 0.447 0.7410], 'Linewidth', 2);
hold on
h2 = plot(TimeIndex(1:per+1),check.BB_liabilities_insolvent_aggregate(1:per+1,1), 'Linewidth', 1);
 
set(get(gca,'YLabel'),'Interpreter','latex','String','Balance Sheets of Insolvent Banks','FontSize',fontsz);
set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',fontsz);
% set(get(gca,'YLabel'),'Rotation',0)
  set(gca, 'Fontsize',13)
hl=legend([h1 h2],{'Assets$$^{k}$$','Liabilities$$^{k}$$'});
 set(hl, 'interpreter', 'latex', 'Location', 'northwest')
 legend('boxoff')
  print -depsc figures/Checks/BB_insolvent;