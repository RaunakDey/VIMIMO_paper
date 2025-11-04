clear all;
clc;

%% datasets.

load('./nophage_control.mat');
load('./triplicate_data.mat');

%% plotting

host1 =host1*1e3;
host2 =host2*1e3;
host3 =host3*1e3;
host4 =host4*1e3;
host5 =host5*1e3;
CBA4 = CBA4*1e3;
CBA18 = CBA18*1e3;
CBA38 = CBA38*1e3;
H100 = H100*1e3;
PSA13_15 = PSA13_15*1e3;

time = time/60;


color_green = [171,193,157]./255;

%% fit growth rates 15.75 hours
g_low = zeros(1,5);
g_high = zeros(1,5);
growthrate = zeros(1,5);

upto_which_point = 28; %used for 6.4 hrs, use 28 for 15.75 hrs.

f = fit(time(1:upto_which_point,:),log(mean(CBA4(1:upto_which_point,:),2)),'poly1');
coefficientValues = coeffvalues(f);
values = confint(cfit(f));
g_low(1) = values(1,1);
g_high(1) = values(2,1);
growthrate(1) = coefficientValues(1);
dev_gr(1) = (g_high(1) - g_low(1))/3; 

f = fit(time(1:upto_which_point,:),log(mean(CBA18(1:upto_which_point,:),2)),'poly1');
coefficientValues = coeffvalues(f);
values = confint(cfit(f));
g_low(2) = values(1,1);
g_high(2) = values(2,1);
growthrate(2) = coefficientValues(1);
dev_gr(2) = (g_high(2) - g_low(2))/3; 


f = fit(time(1:upto_which_point,:),log(mean(CBA38(1:upto_which_point,:),2)),'poly1');
coefficientValues = coeffvalues(f);
values = confint(cfit(f));
g_low(3) = values(1,1);
g_high(3) = values(2,1);
growthrate(3) = coefficientValues(1);
dev_gr(3) = (g_high(3) - g_low(3))/3; 


f = fit(time(1:upto_which_point,:),log(mean(H100(1:upto_which_point,:),2)),'poly1');
coefficientValues = coeffvalues(f);
values = confint(cfit(f));
g_low(4) = values(1,1);
g_high(4) = values(2,1);
growthrate(4) = coefficientValues(1);
dev_gr(4) = (g_high(4) - g_low(4))/3; 


f = fit(time(1:upto_which_point,:),log(mean(PSA13_15(1:upto_which_point,:),2)),'poly1');
coefficientValues = coeffvalues(f);
values = confint(cfit(f));
g_low(5) = values(1,1);
g_high(5) = values(2,1);
growthrate(5) = coefficientValues(1);
dev_gr(5) = (g_high(5) - g_low(5))/3; 



%% fit growth rates 6.4 hours
upto_which_point2 = 10; 

g2_low = zeros(1,5);
g2_high = zeros(1,5);
growthrate2 = zeros(1,5);



f = fit(time(1:upto_which_point2,:),log(mean(CBA4(1:upto_which_point2,:),2)),'poly1');
coefficientValues = coeffvalues(f);
values = confint(cfit(f));
g2_low(1) = values(1,1);
g2_high(1) = values(2,1);
growthrate2(1) = coefficientValues(1);
dev_gr2(1) = (g2_high(1) - g2_low(1))/3; 

f = fit(time(1:upto_which_point2,:),log(mean(CBA18(1:upto_which_point2,:),2)),'poly1');
coefficientValues = coeffvalues(f);
values = confint(cfit(f));
g2_low(2) = values(1,1);
g2_high(2) = values(2,1);
growthrate2(2) = coefficientValues(1);
dev_gr2(2) = (g2_high(2) - g2_low(2))/3; 


f = fit(time(1:upto_which_point2,:),log(mean(CBA38(1:upto_which_point2,:),2)),'poly1');
coefficientValues = coeffvalues(f);
values = confint(cfit(f));
g2_low(3) = values(1,1);
g2_high(3) = values(2,1);
growthrate2(3) = coefficientValues(1);
dev_gr2(3) = (g2_high(3) - g2_low(3))/3; 


f = fit(time(1:upto_which_point2,:),log(mean(H100(1:upto_which_point2,:),2)),'poly1');
coefficientValues = coeffvalues(f);
values = confint(cfit(f));
g2_low(4) = values(1,1);
g2_high(4) = values(2,1);
growthrate2(4) = coefficientValues(1);
dev_gr2(4) = (g2_high(4) - g2_low(4))/3; 


f = fit(time(1:upto_which_point2,:),log(mean(PSA13_15(1:upto_which_point2,:),2)),'poly1');
coefficientValues = coeffvalues(f);
values = confint(cfit(f));
g2_low(5) = values(1,1);
g2_high(5) = values(2,1);
growthrate2(5) = coefficientValues(1);
dev_gr2(5) = (g2_high(5) - g2_low(5))/3; 



%%

load('./../revision/corrected4/v23.mat','chain');
seivd_r = median(chain(40001:50000,28:32));
dev_seivd_r = std(chain(40001:50000,28:32));


     r_onestep = [ 0.18702      0.21341      0.21843      0.27922      0.25593];
     r_onestep_std = [    0.048525     0.066347     0.046714     0.058899     0.052958];

%% plots

figure(2)
subplot(3,3,1)
errorbar(time,mean(host1,2),std(host1'),'LineWidth',2,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor',[1,0,0],Color=[1,0,0]);hold on;
hold on;
errorbar(time,mean(CBA4,2),std(CBA4'), 'LineWidth',2,'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',[0,0,1],Color=[0,0,1] )
set(gca,'FontSize',18);
set(gca,'YScale','log');
ylabel({'Bacterial density';'(cells/mL)'});
xlabel('Time (hrs)');
xticks(0:2:16);
title('CBA 4');

subplot(3,3,2)
errorbar(time,mean(host2,2),std(host2'),'LineWidth',2,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor',[1,0,0],Color=[1,0,0]);hold on;
hold on;
errorbar(time,mean(CBA18,2),std(CBA18'), 'LineWidth',2,'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',[0,0,1],Color=[0,0,1] )
set(gca,'FontSize',18);
set(gca,'YScale','log');
xlabel('Time (hrs)');
xticks(0:2:16);
title('CBA 18');



subplot(3,3,3)
errorbar(time,mean(host3,2),std(host3'),'LineWidth',2,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor',[1,0,0],Color=[1,0,0]);hold on;
hold on;
errorbar(time,mean(CBA38,2),std(CBA38'), 'LineWidth',2,'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',[0,0,1],Color=[0,0,1] )
set(gca,'FontSize',18);
set(gca,'YScale','log');
xlabel('Time (hrs)');
xticks(0:2:16);
title('CBA 38');



subplot(3,3,4)
errorbar(time,mean(host4,2),std(host4'),'LineWidth',2,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor',[1,0,0],Color=[1,0,0]);hold on;
hold on;
errorbar(time,mean(H100,2),std(H100'), 'LineWidth',2,'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',[0,0,1],Color=[0,0,1] )
set(gca,'FontSize',18);
set(gca,'YScale','log');
ylabel({'Bacterial density';'(cells/mL)'});
xlabel('Time (hrs)');
xticks(0:2:16);
title('PSA H100');


subplot(3,3,5)
errorbar(time,mean(host5,2),std(host5'),'LineWidth',2,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor',[1,0,0],Color=[1,0,0]);hold on;
hold on;
errorbar(time,mean(PSA13_15,2),std(PSA13_15'), 'LineWidth',2,'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',[0,0,1],Color=[0,0,1] )
set(gca,'FontSize',18);
set(gca,'YScale','log');
xlabel('Time (hrs)');
xticks(0:2:16);
title('PSA 13-15');


subplot(3,3,6)
plot(time,mean(host1+host2+host3+host4+host5,2),'LineWidth',2,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor',[1,0,0],Color=[1,0,0]);hold on;
hold on;
plot(time,mean(CBA4+CBA18+CBA38+ H100+PSA13_15,2), 'LineWidth',2,'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',[0,0,1],Color=[0,0,1] )
set(gca,'FontSize',18);
set(gca,'YScale','log');
xlabel('Time (hrs)');
xticks(0:2:16);
title('All hosts');


subplot(3,3,[7,8,9])

errorbar((1:5)-0.15,growthrate2,dev_gr2,'Marker','square',MarkerSize = 8, MarkerFaceColor='k',Color='k',LineStyle='none',LineWidth=2);
hold on;
errorbar((1:5)-0.05,growthrate,dev_gr,'Marker','square',MarkerSize = 8, MarkerFaceColor='b',Color='b',LineStyle='none',LineWidth=2);
errorbar((1:5)+0.05,seivd_r,dev_seivd_r,'Marker','^',MarkerSize = 8, MarkerFaceColor='r',Color='r',LineStyle='none',LineWidth=2);
xticks(1:5);
errorbar((1:5)+0.15,[0.19, 0.245,0.22,0.28,0.25],[0.019,0.021,0.030,0.052,0.044],'Marker','^',MarkerSize = 8,MarkerFaceColor=[255, 219, 88]./255, MarkerEdgeColor = [255, 219, 88]./255, Color =[255, 219, 88]./255,  LineStyle='none',LineWidth=2)
%errorbar((1:5)+0.15, r_onestep,r_onestep_std, 'o','MarkerSize',8,'MarkerFaceColor',color_green,'Color',color_green,'LineWidth',2)
xticklabels({'CBA4','CBA18','CBA38','PSA H100','PSA 13-15'})
set(gca,'FontSize',18);
legend('No phage control qPCR fitted 5.25 hours','No phage control qPCR fitted to the complete 15.75 hours','With phages qPCR (fitted with SEIVD model) to 5\times5 community','Doubling time OD experiments','Bayesian estimate from pairwise experiments');
ylabel({'Bacterial growth rates' ; '(cells/hr)'})
xlim([0.5,5.5]);
grid on;
