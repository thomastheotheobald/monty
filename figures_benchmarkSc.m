%% Load the Benchmark Scenario and create the figures
load Workspaces/BenchmarkWS/21112017/MC2_workspace1000.mat;

% aktueller Workspace: 'Workspaces/WS/WS02082017_1000'
%% Load the confidence bands
% It is neccessary to run the MonteCarlo file first
load ConfidenceBands.mat
%% 
clc;
close all;

%%
% Declare fontsize
fsize = 15;


%% Histograms of SPV bankruptcies

norm    = 10;
xlimit    = 100;
ylimit = 0.135;
xlimhist = 0.09;


% % Selected Figures of the Benchmark Case
%% The equity Ratios
c = 30;
a = 0.0;
b = 0.7;
j_ind = a:(b-a)/c:b;
j_ind = j_ind';


% Time horizon of the plots
 per = 51;
 
 % Choose the representative replication
 rep = 9; 
%% Plot the equity ratios and the SPV profits into one figure
figure('Visible','off')
XX = 1:100; XX = XX';
for j = 1:c
    if j ==1
  [ax,h1,h2] = plotyy(XX,MC.ER_k(:,j,rep), XX, MC.pi_spv(:,rep))  ;
  set(ax(1),{'ycolor'},{'k'},{'YLim'},{[0.09,0.135]},{'YTick'},{0.09:0.01:0.135})
  set(ax(2),{'ycolor'},{'k'},{'YLim'},{[-2*(10^(-3))*2+0.001 1*(10^(-3))*2+0.001]},'FontSize',13)
  set(h1,'color',[0.6 0.6 0.6], 'Linewidth', 2) 

  set(h2,'color',[0 0.447 0.7410], 'Linewidth', 2) 

    else
  plot(MC.ER_k(:,j,rep),'color', [j_ind(j) j_ind(j) j_ind(j)], 'Linewidth', 1.5)
    end
  hold on
end
  plot([0,200],[ER_T, ER_T],':k','Linewidth', 1);   
  set(gca, 'Fontsize',13)
  xlim([0,xlimit]);
  ylim([0.09,0.135])
    set(get(gca,'YLabel'),'Interpreter','latex','String','$$ER^{k}_{t}$$'...
        ,'Units', 'Normalized', 'Position', [-0.1, 0.46, 0],'FontSize',18);
    set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',18);
    % Rotate Y-Label
    set(get(gca,'YLabel'),'Rotation',0)
  print -depsc figures/Benchmark/ER;

% SPV profits
  figure('Visible','off')
  plot(MC.pi_spv(:,rep), 'Linewidth', 1.5)
  hold on
  plot([0,200],[0, 0],':k','Linewidth', 1);   
  set(gca, 'Fontsize',13)
  xlim([0,xlimit]);
%   ylim([0.09,ylimit]);
  print -depsc figures/Benchmark/pi_spv;  
  
  
%% The economic growth rate
  TimeIndex = 1:N;
  TimeIndex = (TimeIndex-1)';
  
  
  figure('Visible','off')
X = res.g_Cup-res.g_median + MC.g(:,rep);
Y = MC.g(:,rep)-(res.g_median-res.g_Cbot);
%   plot and patch the corresponding confidence intervals
fill([TimeIndex(1:per)' fliplr(TimeIndex(1:per)')],[X(1:per)' fliplr(Y(1:per)')],[0.9 0.9 0.9])
hold on
plot(TimeIndex(1:per),MC.g(1:per,rep),'color',[0 0.447 0.7410], 'Linewidth', 1.5)
  ylim([0.015,0.05]);
set(get(gca,'YLabel'),'Interpreter','latex','String','$$g_{t}$$'...
        ,'Units', 'Normalized', 'Position', [-0.1, 0.46, 0],'FontSize',18); 
set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',18);
set(get(gca,'YLabel'),'Rotation',0)
  set(gca, 'Fontsize',11)
  print -depsc figures/Benchmark/g;

  %% The associated output evolution 
   figure('Visible','off')
   X = res.y_Cup-res.y_median + MC.y(:,rep);
   Y = MC.y(:,rep)-(res.y_median-res.y_Cbot);
%   plot and patch the corresponding confidence intervals
fill([TimeIndex(1:per)' fliplr(TimeIndex(1:per)')],[X(1:per)' fliplr(Y(1:per)')],[0.9 0.9 0.9])
hold on
  plot(TimeIndex(1:per),MC.y(1:per,rep),'color',[0 0.447 0.7410], 'Linewidth', 1.5)
set(get(gca,'YLabel'),'Interpreter','latex','String','$$y_{t}$$'...
        ,'Units', 'Normalized', 'Position', [-0.1, 0.46, 0],'FontSize',18); 
    set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',18);
        ylim([50,250])
    % Rotate Y-Label
set(get(gca,'YLabel'),'Rotation',0)
  set(gca, 'Fontsize',11)
  print -depsc figures/Benchmark/y; 
 %% Plot both in one figure with subfigures
 
  figure('Visible','off')
  subplot(2,1,1)
X = res.y_Cup-res.y_median + MC.y(:,rep);
Y = MC.y(:,rep)-(res.y_median-res.y_Cbot);
%   plot and patch the corresponding confidence intervals
fill([TimeIndex(1:per)' fliplr(TimeIndex(1:per)')],[X(1:per)' fliplr(Y(1:per)')],[0.9 0.9 0.9])
hold on
plot(TimeIndex(1:per),MC.y(1:per,rep),'color',[0 0.447 0.7410], 'Linewidth', 1.5)
set(get(gca,'YLabel'),'Interpreter','latex','String','$$y_{t}$$'...
        ,'Units', 'Normalized', 'Position', [-0.1, 0.46, 0],'FontSize',18); 
    set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',18);
        ylim([50,250])
    % Rotate Y-Label
set(get(gca,'YLabel'),'Rotation',0)
  set(gca, 'Fontsize',11)

subplot(2,1,2)
X = res.g_Cup-res.g_median + MC.g(:,rep);
Y = MC.g(:,rep)-(res.g_median-res.g_Cbot);
%   plot and patch the corresponding confidence intervals
fill([TimeIndex(1:per)' fliplr(TimeIndex(1:per)')],[X(1:per)' fliplr(Y(1:per)')],[0.9 0.9 0.9])
hold on
plot(TimeIndex(1:per),MC.g(1:per,rep),'color',[0 0.447 0.7410], 'Linewidth', 1.5)
  ylim([0.015,0.05]);
set(get(gca,'YLabel'),'Interpreter','latex','String','$$g_{t}$$'...
        ,'Units', 'Normalized', 'Position', [-0.1, 0.46, 0],'FontSize',18); 
    set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',18);
    % Rotate Y-Label
set(get(gca,'YLabel'),'Rotation',0)
  set(gca, 'Fontsize',11)
  print -depsc figures/Benchmark/yg; 
  

%% Bankrupt firms and banks

figure('Visible','off')
X = res.bankrupt_banks_Cup-res.bankrupt_banks_median + MC.bankrupt_banks(:,rep);
Y = MC.bankrupt_banks(:,rep)-(res.bankrupt_banks_median-res.bankrupt_banks_Cbot);
%   plot and patch the corresponding confidence intervals
fill([TimeIndex(1:per)' fliplr(TimeIndex(1:per)')],[X(1:per)' fliplr(Y(1:per)')],[0.9 0.9 0.9])
hold on
plot(TimeIndex(1:per),MC.bankrupt_banks(1:per,rep),'color',[0 0.447 0.7410], 'Linewidth', 1.5)
  xlim([0,per-1]);
set(get(gca,'YLabel'),'Interpreter','latex','String','insolvent banks (in \% of $$N^{B}$$)','FontSize',18); 
set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',18);
  set(gca, 'Fontsize',11)
  print -depsc figures/Benchmark/bankrupt_banks;

figure('Visible','off')
X = res.bankrupt_firms_Cup-res.bankrupt_firms_median + MC.bankrupt_firms(:,rep);
Y = MC.bankrupt_firms(:,rep)-(res.bankrupt_firms_median-res.bankrupt_firms_Cbot);
%   plot and patch the corresponding confidence intervals
fill([TimeIndex(1:per)' fliplr(TimeIndex(1:per)')],[X(1:per)' fliplr(Y(1:per)')],[0.9 0.9 0.9])
hold on
plot(TimeIndex(1:per),MC.bankrupt_firms(1:per,rep),'color',[0 0.447 0.7410], 'Linewidth', 1.5)
%   ylim([0.015,0.05]);
set(get(gca,'YLabel'),'Interpreter','latex','String','insolvent firms (in \% of $$N^{C}$$)','FontSize',18); 
    set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',18);
  set(gca, 'Fontsize',11)
  print -depsc figures/Benchmark/bankrupt_firms;
  
%%  Putting all together into one matrix of plots

  figure('Visible','off')
    subplot(2,2,1)
X = res.bankrupt_firms_Cup-res.bankrupt_firms_median + MC.bankrupt_firms(:,rep);
Y = MC.bankrupt_firms(:,rep)-(res.bankrupt_firms_median-res.bankrupt_firms_Cbot);
%   plot and patch the corresponding confidence intervals
fill([TimeIndex(1:per)' fliplr(TimeIndex(1:per)')],[X(1:per)' fliplr(Y(1:per)')],[0.9 0.9 0.9])
hold on
plot(TimeIndex(1:per),MC.bankrupt_firms(1:per,rep),'color',[0 0.447 0.7410], 'Linewidth', 1.5)
%   ylim([0.015,0.05]);
set(get(gca,'YLabel'),'Interpreter','latex','String','insolvent firms (in \% of $$N^{C}$$)','FontSize',18); 
    set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',18);
      xlim([0,per-1]);
  set(gca, 'Fontsize',11)
  
     subplot(2,2,2)
X = res.bankrupt_banks_Cup-res.bankrupt_banks_median + MC.bankrupt_banks(:,rep);
Y = MC.bankrupt_banks(:,rep)-(res.bankrupt_banks_median-res.bankrupt_banks_Cbot);
%   plot and patch the corresponding confidence intervals
fill([TimeIndex(1:per)' fliplr(TimeIndex(1:per)')],[X(1:per)' fliplr(Y(1:per)')],[0.9 0.9 0.9])
hold on
plot(TimeIndex(1:per),MC.bankrupt_banks(1:per,rep),'color',[0 0.447 0.7410], 'Linewidth', 1.5)
%   ylim([0.015,0.05]);
set(get(gca,'YLabel'),'Interpreter','latex','String','insolvent banks (in \% of $$N^{B}$$)','FontSize',18); 
set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',18);
  set(gca, 'Fontsize',11)
    xlim([0,per-1]);
  
  subplot(2,2,3)
X = res.y_Cup-res.y_median + MC.y(:,rep);
Y = MC.y(:,rep)-(res.y_median-res.y_Cbot);
%   plot and patch the corresponding confidence intervals
fill([TimeIndex(1:per)' fliplr(TimeIndex(1:per)')],[X(1:per)' fliplr(Y(1:per)')],[0.9 0.9 0.9])
hold on
plot(TimeIndex(1:per),MC.y(1:per,rep),'color',[0 0.447 0.7410], 'Linewidth', 1.5)
set(get(gca,'YLabel'),'Interpreter','latex','String','$$y_{t}$$'...
        ,'Units', 'Normalized', 'Position', [-0.1, 0.46, 0],'FontSize',18); 
    set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',18);
    % Rotate Y-Label
    ylim([50,250])
set(get(gca,'YLabel'),'Rotation',0)
  set(gca, 'Fontsize',11)
ylab = get(gca,'YLabel'); 
set(ylab,'Position',get(ylab,'Position') + [0 .1 0])
  xlim([0,per-1]);

subplot(2,2,4)
X = res.g_Cup-res.g_median + MC.g(:,rep);
Y = MC.g(:,rep)-(res.g_median-res.g_Cbot);
%   plot and patch the corresponding confidence intervals
fill([TimeIndex(1:per)' fliplr(TimeIndex(1:per)')],[X(1:per)' fliplr(Y(1:per)')],[0.9 0.9 0.9])
hold on
plot(TimeIndex(1:per),MC.g(1:per,rep),'color',[0 0.447 0.7410], 'Linewidth', 1.5)
  ylim([0.012,0.05]);
set(get(gca,'YLabel'),'Interpreter','latex','String','$$g_{t}$$'...
        ,'Units', 'Normalized', 'Position', [-0.1, 0.46, 0],'FontSize',18); 
    set(get(gca,'XLabel'),'Interpreter','latex','String','time','FontSize',18);
    % Rotate Y-Label
set(get(gca,'YLabel'),'Rotation',0)
  set(gca, 'Fontsize',11)
    xlim([0,per-1]);
    ylab = get(gca,'YLabel'); 
set(ylab,'Position',get(ylab,'Position') + [0 .09 0])
  print -depsc figures/Benchmark/4panels; 