%% Illustration

%% Define the framework
figure(1)

subplot(3,3,1)
plot(ma.m_d)
axis([15 30 50 250])
title('deposits (level)','FontSize',10)

subplot(3,3,2)
plot(ma.pi_spvA)
axis([15 30 -0.1 0.1])
title('spv profits','FontSize',10)

subplot(3,3,3)
plot(MC.bankrupt_banks)
axis([15 30 0 25])
title('bank insolvencies (%)','FontSize',10) 

subplot(3,3,4)
plot(ma.l_c)
axis([15 30 50 250])
title('credit (level)','FontSize',10)

subplot(3,3,5)
plot(ma.K_solvent)
axis([15 30 150 600])
title('capital stock (level)','FontSize',10)

subplot(3,3,6)
plot(ma.i)
axis([15 30 0 15])
title('net investment (level)','FontSize',10)

subplot(3,3,7)
plot(ma.c)
axis([15 30 30 70])
title('consumption(level)','FontSize',10)

subplot(3,3,8)
plot(ma.y)
axis([15 30 30 80])
title('GDP (level)','FontSize',10)

subplot(3,3,9)
plot(ma.growth)
axis([15 30 -0.06 0.1])
title('GDP growth (%)','FontSize',10)
