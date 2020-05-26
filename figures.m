%% Illustration

%% Growth Rate


% [div.g_sort, div.g_order] = sort(MC.g_std);
% fig.g_min = MC.g(:,div.g_order(:,1));
% fig.g_max = MC.g(:,div.g_order(:,end));
% fig.g     = MC.g(:,div.g_order(:,Ncol/2));
% 
% figure
% plot(fig.g,'Linewidth',2)
% hold on
% plot(fig.g_min,'--','Linewidth',1)
% hold on
% plot(fig.g_max,'--','Linewidth',1)

%%
% Declare fontsize
fsize = 15;

% % Monte Carlo simulation results

% figure('Visible','off')
% plot(res.g_mean*100,'b','Linewidth',2)
% hold on
% plot(res.g_Cup*100,'--k','Linewidth',1)
% hold on
% plot(res.g_Cbot*100,'--k','Linewidth',1)
% set(get(gca,'XLabel'),'Interpreter','latex','String','$$time$$','FontSize',fsize); 
% set(get(gca,'YLabel'),'Interpreter','latex','String','$$g(\%)$$','FontSize',fsize,'Rotation',0); 
% set(gca, 'Fontsize',11);
%  ylab = get(gca,'YLabel'); 
% set(ylab,'Position',get(ylab,'Position') - [0 .04 0])
% box off

%% Histograms of SPV bankruptcies

norm    = 10;
xlimit    = 100;
ylimit = 0.135;
xlimhist = 0.09;
% % TimeIndex = 1:xlim+20;
% %     TimeIndex = TimeIndex -1;
% Vor SPV Pleite
if replications >1
        if q == 1 %Requires 'Results.m'
            figure('Visible','off')
                hist(MC.ER_k(MC.spv_bankrupt -1,:));
                xlim([xlimhist,0.17]);
                set(gca, 'Fontsize',13)
                print -depsc figures/MC1/histER1;
                
            figure('Visible','off')
                hist(MC.ER_k(MC.spv_bankrupt,:));
                xlim([xlimhist,0.17]);
                set(gca, 'Fontsize',13)
                print -depsc figures/MC1/histER2;
                
            figure('Visible','off')
                plot(MC.ER_k(:,:,1), 'Linewidth', 2)
                hold on
                plot([0,200],[ER_T, ER_T],':k','Linewidth', 1);                
                set(gca, 'Fontsize',13)
                xlim([0,xlimit]);
                print -depsc figures/MC1/ER;

        elseif q == 2
            figure('Visible','off')
                hist(MC.ER_k(MC.spv_bankrupt -1,:));
                xlim([xlimhist,0.17]);
                set(gca, 'Fontsize',13)
            print -depsc figures/MC2/histER1;
            
            figure('Visible','off')
                hist(MC.ER_k(MC.spv_bankrupt,:));
                xlim([xlimhist,0.17]);
                set(gca, 'Fontsize',13)
                print -depsc figures/MC2/histER2;

            figure('Visible','off')
                plot(MC.ER_k(:,:,1), 'Linewidth', 2)
                hold on
                plot([0,200],[ER_T, ER_T],':k','Linewidth', 1); 
                xlim([0,xlimit]);
                set(gca, 'Fontsize',13)
                print -depsc figures/MC2/ER;   
                
            figure('Visible','off')
                plot(MC.ER_k(:,:,1), 'Linewidth', 2)
                hold on
                plot(MC.spvA(:,1)/norm, '--k','Linewidth', 2)
                hold on
                plot([0,200],[ER_T, ER_T],':k','Linewidth', 1);
                xlim([0,xlimit])                
                set(gca, 'Fontsize',13)
                print -depsc figures/MC2/ER_spvA; 
        elseif q == 3
            figure('Visible','off')
                hist(MC.ER_k(MC.spv_bankrupt -1,:));
                xlim([xlimhist,0.17]);
                set(gca, 'Fontsize',13)
            print -depsc figures/MC3/histER1;
            
            figure('Visible','off')
                hist(MC.ER_k(MC.spv_bankrupt,:));
                xlim([xlimhist,0.17]);
                set(gca, 'Fontsize',13)
                print -depsc figures/MC3/histER2;
                
            figure('Visible','off')
                plot(MC.ER_k(:,:,1), 'Linewidth', 2)
                hold on
                plot([0,200],[ER_T, ER_T],':k','Linewidth', 1); 
                xlim([0,xlimit]);                
                set(gca, 'Fontsize',13)
                print -depsc figures/MC3/ER; 
                
            figure('Visible','off')
                plot(MC.ER_k(:,:,1), 'Linewidth', 2)
                hold on
                plot(MC.spvA(:,1)/norm, '--k','Linewidth', 2)
                hold on
                plot([0,200],[ER_T, ER_T],':k','Linewidth', 1);
                xlim([0,xlimit])                
                set(gca, 'Fontsize',13)
                print -depsc figures/MC3/ER_spvA; 
                
            
% % %                 %% Transaction costs
% % %             figure('Visible','off')
% % %                 plot(MC.transaction_average_single(:,1,1), 'k','Linewidth',1)
% % %                 hold on
% % %                 plot(MC.transaction_average_single(:,1,2),'--k' ,'Linewidth', 2)
% % %                 hold on
% % %                 plot(MC.transaction_average_single(:,1,3),'color', [0.5 0.5 0.5],'Linewidth', 2)
% % %                 set(gca, 'Fontsize',13)
% % %                 print -depsc figures/transactioncosts;
% % %                 
% % %             figure('Visible','off')
% % %                 plot(MC.transaction_average_all(:,1), 'k','Linewidth',1)
% % %                 hold on
% % %                 plot(MC.transaction_average_all(:,2),'--k' ,'Linewidth', 2)
% % %                 hold on
% % %                 plot(MC.transaction_average_all(:,3),'color', [0.5 0.5 0.5],'Linewidth', 2)
% % %                 set(gca, 'Fontsize',13)
% % %                 print -depsc figures/transactioncosts_all;
% % % 
% % % % Using the mean                
% % %              figure('Visible','off')
% % %                 plot(MC.transaction_average_all(:,1), 'k','Linewidth',1)
% % %                 hold on
% % %                 plot(MC.transaction_average_all(:,2),'color', [0.2 0.2 0.2],'Linewidth', 2)
% % %                 hold on
% % %                 plot(MC.transaction_average_all(:,3),'color', [0.5 0.5 0.5],'Linewidth', 2)
% % %                 hold on
% % %                 plot(MC.bankrupt_bank_average(:,1), ':k','Linewidth',1)
% % %                 hold on
% % %                 plot(MC.bankrupt_bank_average(:,2),':' ,'color', [0.2 0.2 0.2],'Linewidth',2)
% % %                 hold on
% % %                 plot(MC.bankrupt_bank_average(:,3),':','color', [0.5 0.5 0.5],'Linewidth',2)
% % %                 set(gca, 'Fontsize',13)
% % %                 print -depsc figures/transactioncosts_bankruptcies;
% % %                 
% % % % Using the median
% % %              figure('Visible','off')
% % %                 plot(MC.transaction_median_all(:,1), 'k','Linewidth',1)
% % %                 hold on
% % %                 plot(MC.transaction_median_all(:,2),'color', [0.2 0.2 0.2],'Linewidth', 2)
% % %                 hold on
% % %                 plot(MC.transaction_median_all(:,3),'color', [0.5 0.5 0.5],'Linewidth', 2)
% % %                 hold on
% % %                 plot(MC.bankrupt_bank_median(:,1), ':k','Linewidth',1)
% % %                 hold on
% % %                 plot(MC.bankrupt_bank_median(:,2),':' ,'color', [0.2 0.2 0.2],'Linewidth',2)
% % %                 hold on
% % %                 plot(MC.bankrupt_bank_median(:,3),':','color', [0.5 0.5 0.5],'Linewidth',2)
% % %                 set(gca, 'Fontsize',13)
% % %                 print -depsc figures/transactioncosts_median_bankruptcies;
        end
end