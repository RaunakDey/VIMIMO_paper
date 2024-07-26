
load('./../revision/corrected4/23b.mat');
close all;

D2 = D_after;
Dc = pars_afterinf.Dc ;
Dc2 = pars_afterinf.Dc2;
Dc3 = pars_afterinf.Dc3 ;
Dc4 = pars_afterinf.Dc4;
Dc5= pars_afterinf.Dc5;


red1 = [1 0.1 0.1];
red2 = [0.8 0.3 0.3];
red3 = [0.6 0.4 0.4];

blue1 = [0.1 0.1 1];
blue2 = [0.3 0.3 0.8];
%%
figure(1)
subplot(1,2,1)





prob = 1 - 1./(1+(D2./Dc).^2);
prob2 = 1 - 1./(1+(D2./Dc2).^2);
prob3 = 1 - 1./(1+(D2./Dc3).^2);
prob4 = 1 - 1./(1+(D2./Dc4).^2);
prob5 = 1 - 1./(1+(D2./Dc5).^2); 


plot(t2,prob,'LineWidth',2.5,'Color',red1);
hold on; 
plot(t2,prob2,'LineWidth',2.5,'Color',red2);
plot(t2,prob3,'LineWidth',2.5,'Color',red3);
plot(t2,prob4,'LineWidth',2.5,'Color',blue1);
plot(t2,prob5,'LineWidth',2.5,'Color',blue2);
xlabel('Time (hours)');
ylabel('Probability of lysis inhibition');
set(gca,'FontSize',24);
xlim([0 8]);
xticks(0:2:16)
legend('CBA 4', 'CBA 18', 'CBA 38', 'PSA H100', 'PSA 13-15');

subplot(1,2,2)
plot(D2,prob,'LineWidth',2.5,'Color',red1);
hold on;
plot(D2,prob2,'LineWidth',2.5,'Color',red2);
plot(D2,prob3,'LineWidth',2.5,'Color',red3);
plot(D2,prob4,'LineWidth',2.5,'Color',blue1);
plot(D2,prob5,'LineWidth',2.5,'Color',blue2);

set(gca, 'YScale', 'log')
xlabel({'Debris concentration', '(10^7 Lysed cells/ml)'});
ylabel({'Probability of lysis inhibition'});
set(gca,'FontSize',24);
yline(0.5,'--k',LineWidth=1.5);

xline(Dc, '--', LineWidth= 1,Color=red1);
xline(Dc2, '--', LineWidth= 1,Color=red2);
xline(Dc3, '--', LineWidth= 1,Color=red3);
xline(Dc4, '--', LineWidth= 1,Color=blue1);
xline(Dc5, '--', LineWidth= 1,Color=blue2);


%xline(Dc,'--r',LineWidth=1.5);
%yline(0.5,'--b',LineWidth=1.5);
%xline(Dc2,'--b',LineWidth=1.5);
xlim([0 6e7])


xticks([0,2e7,4e7,6e7,8e7,10e7]);
xticklabels({'0','2','4','6','8','10'});
yticks([0.001,0.01,0.1,0.5,1]);
ylim([0.001 1])
legend('CBA 4', 'CBA 18', 'CBA 38', 'PSA H100', 'PSA 13-15');


% subplot(1,3,3)
% set(gca, 'YScale', 'linear')
% plot(t2,D2,'LineWidth',2.5,'Color','r');
% xlabel('Time (hours)');
% xlim([0 16]);
% xticks(0:2:16)
% ylabel({'Debris concentration', '(10^7 Lysed cells/ml)'});
% set(gca,'FontSize',24);