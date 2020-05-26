

div.sortx      = zeros(N,Ncol);
div.orderx     = zeros(N,Ncol);

res.g_Cup       = zeros(N,1);
res.g_Cbot      = zeros(N,1);


div.Cup = 0.84;
div.Cbot = 0.16;

div.Cup_index = ceil(div.Cup*Ncol);
div.Cbot_index = ceil(div.Cbot*Ncol);

   res.g_mean =  mean(MC.g,2);
   res.g_median = nanmedian(MC.g,2);
   res.y_median = nanmedian(MC.y,2);
   res.bankrupt_banks_median = nanmedian(MC.bankrupt_banks,2);
   res.bankrupt_firms_median = nanmedian(MC.bankrupt_firms,2);
for s =1 : N   
   [div.sortx(s,:),div.orderx(s,:)] = sort(MC.g(s,:));
   res.g_Cup(s,:)  = div.sortx(s,div.Cup_index); 
   res.g_Cbot(s,:)  = div.sortx(s,div.Cbot_index);
   
   [div.sortx(s,:),div.orderx(s,:)] = sort(MC.y(s,:));
   res.y_Cup(s,:)  = div.sortx(s,div.Cup_index); 
   res.y_Cbot(s,:)  = div.sortx(s,div.Cbot_index);
   
   [div.sortx(s,:),div.orderx(s,:)] = sort(MC.bankrupt_banks(s,:));
   res.bankrupt_banks_Cup(s,:)  = div.sortx(s,div.Cup_index); 
   res.bankrupt_banks_Cbot(s,:)  = div.sortx(s,div.Cbot_index); 
   
   [div.sortx(s,:),div.orderx(s,:)] = sort(MC.bankrupt_firms(s,:));
   res.bankrupt_firms_Cup(s,:)  = div.sortx(s,div.Cup_index); 
   res.bankrupt_firms_Cbot(s,:)  = div.sortx(s,div.Cbot_index); 
end

save ConfidenceBands res;
