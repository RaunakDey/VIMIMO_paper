clear all;
clc;


%% control

cba4_nophage =  1e3*[107.128	154.3845	94.866;
2849.473333	3863.433333	4287.333333];


cba18_nophage = 1e3*[235.3305	243.9835	348.4855;
87981.66667	18953.03333	3613.233333];


cba38_nophage  = 1e3*[122.2875	111.594	135.1055;
12648.2	15114.36667	13963.13333];

h100_nophage = 1e3*[94.441	5.05	68.0175;
2678.223333	20936.93333	5679];

psa1315_nophage = 1e3*[354.129	102.573	296.4605
26023.2	24411.66667	48637.66667];


%% 1
cba18_2_18_host = 1e3*[4.30E+02	2.66E+02	4.71E+02;
2.90E+03	2.27E+03	2.71E+03];

cba18_2_18_phage = 1e3*[1.50E+04	1.50E+04	1.37E+04;
6.01E+08	5.57E+08	5.58E+08];

X ={'Uninfected host','Infected host','Phage'};

figure(2)
subplot(3,3,1)
what_to_plot = [mean(cba18_nophage');mean(cba18_2_18_host'); mean(cba18_2_18_phage')];
bar(X,what_to_plot);
set(gca,'YScale','log');
ylim([1e3 1e12]);
ylabel('Genome copies/ml');
set(gca,'FontSize',20);
yticks([1e3 1e5 1e7 1e9 1e11]);
title('CBA 18:2 on CBA 18');


%% 2
cba18_3_4_host = 1e3*[1.04E+02	1.14E+02	2.28E+02;
5.37E+02	5.26E+02	5.16E+02];

cba18_3_4_phage = 1e3*[1.11E+04	1.30E+04	1.29E+04;
1.97E+07	1.91E+07	1.90E+07];

subplot(3,3,2)
what_to_plot = [mean(cba4_nophage');mean(cba18_3_4_host'); mean(cba18_3_4_phage')];
bar(X,what_to_plot);
set(gca,'YScale','log');
ylim([1e3 1e12]);
ylabel('Genome copies/ml');
set(gca,'FontSize',20);
yticks([1e3 1e5 1e7 1e9 1e11]);
title('CBA 18:3 on CBA 4');

%% 3
cba18_3_18_host = 1e3*[3.32E+02	2.26E+02	3.80E+02;
1.45E+02	1.02E+02	7.51E+01];

cba18_3_18_phage = 1e3*[1.08E+04	1.12E+04	1.22E+04;
6.36E+07	5.93E+07	4.98E+07];


subplot(3,3,3)
what_to_plot = [mean(cba18_nophage');mean(cba18_3_18_host'); mean(cba18_3_18_phage')];
bar(X,what_to_plot);
set(gca,'YScale','log');
ylim([1e3 1e12]);
ylabel('Genome copies/ml');
set(gca,'FontSize',20);
yticks([1e3 1e5 1e7 1e9 1e11]);
title('CBA 18:3 on CBA 18');

%% 4

cba38_1_18_host = 1e3*[2.35E+02	3.26E+02	3.33E+02;
3.87E+03	3.29E+03	3.80E+03];

cba38_1_18_phage = 1e3*[4.32E+05	4.75E+05	4.76E+05; %this had a mistake
4.48E+05	4.46E+05	4.00E+05];

subplot(3,3,4)
what_to_plot = [mean(cba18_nophage');mean(cba38_1_18_host'); mean(cba38_1_18_phage')];
bar(X,what_to_plot);
set(gca,'YScale','log');
ylim([1e3 1e12]);
ylabel('Genome copies/ml');
set(gca,'FontSize',20);
yticks([1e3 1e5 1e7 1e9 1e11]);
title('CBA 38:1 on CBA 18');




%% 5 

cba38_1_38_host = 1e3*[1.07E+02	1.15E+02	1.00E+02;
1.28E+02	8.30E+01	2.68E+02];

cba38_1_38_phage = 1e3*[3.86E+05	4.23E+05	4.22E+05;
9.79E+06	9.57E+06	9.60E+06];


subplot(3,3,5)
what_to_plot = [mean(cba38_nophage');mean(cba38_1_38_host'); mean(cba38_1_38_phage')];
bar(X,what_to_plot);
set(gca,'YScale','log');
ylim([1e3 1e12]);
ylabel('Genome copies/ml');
set(gca,'FontSize',20);
yticks([1e3 1e5 1e7 1e9 1e11]);
title('CBA 38:1 on CBA 38');

%% 6

hp1_h100_host = 1e3*[8.79E+01	5.59E+01	2.70E+01;
3.05E+01	3.60E+00	2.01E+01];

hp1_h100_phage = 1e3*[2.55E+03	3.41E+03	3.51E+03;
8.25E+07	7.98E+07	7.45E+07];

subplot(3,3,6)
what_to_plot = [mean(h100_nophage');mean(hp1_h100_host'); mean(hp1_h100_phage')];
bar(X,what_to_plot);
set(gca,'YScale','log');
ylim([1e3 1e12]);
ylabel('Genome copies/ml');
set(gca,'FontSize',20);
yticks([1e3 1e5 1e7 1e9 1e11]);
title('PSA HP1 on PSA H100');

%% 7
hp1_1315_host = 1e3*[1.58E+02	9.28E+01	2.59E+02;
8.97E+00	7.95E+00	1.01E+01];

hp1_1315_phage = 1e3*[4.22E+03	3.60E+03	3.37E+03;
5.54E+07	5.31E+07	5.46E+07];

subplot(3,3,7)
what_to_plot = [mean(psa1315_nophage');mean(hp1_1315_host'); mean(hp1_1315_phage')];
bar(X,what_to_plot);
set(gca,'YScale','log');
ylim([1e3 1e12]);
ylabel('Genome copies/ml');
set(gca,'FontSize',20);
yticks([1e3 1e5 1e7 1e9 1e11]);
title('PSA HP1 on PSA 13-15');

%% 8

hs6_h100_host = 1e3*[1.16E+01	6.26E+01	6.66E+01;
9.89E+02	2.36E+02	9.40E+01];

hs6_h100_phage = 1e3*[1.05E+04	6.57E+06	1.09E+04;
2.08E+08	1.68E+08	1.85E+08];

subplot(3,3,8)
what_to_plot = [mean(h100_nophage');mean(hs6_h100_host'); mean(hs6_h100_phage')];
bar(X,what_to_plot);
set(gca,'YScale','log');
ylim([1e3 1e12]);
ylabel('Genome copies/ml');
set(gca,'FontSize',20);
yticks([1e3 1e5 1e7 1e9 1e11]);
title('PSA HS6 on PSA H100');

%% 9
hs6_1315_host = 1e3*[9.02E+01	4.22E+02	1.03E+02;
6.25E+02	7.06E+02	1.30E+03];

hs6_1315_phage = 1e3*[1.02E+04	1.07E+04	1.21E+04;
1.03E+08	9.02E+07	9.36E+07];

subplot(3,3,9)
what_to_plot = [mean(psa1315_nophage');mean(hs6_1315_host'); mean(hs6_1315_phage')];
bar(X,what_to_plot);
set(gca,'YScale','log');
ylim([1e3 1e12]);
ylabel('Genome copies/ml');
set(gca,'FontSize',20);
yticks([1e3 1e5 1e7 1e9 1e11]);
title('PSA HS6 on PSA 13-15');

% %%
% figure()
% 
% subplot(5,4,1)
% errorbar(time,mean(cba18_2_18_phage'),std(cba18_2_18_phage'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e12]);
% yticks([1e1 1e3 1e5 1e7 1e9]);
% xlim([-1 17]);
% title('cba18 2 on 18 phage')
% ylabel('genome copies/ml');
% 
% 
% 
% subplot(5,4,2)
% errorbar(time,mean(cba18_2_18_host'),std(cba18_2_18_host'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e8]);
% yticks([1e1 1e2 1e3 1e4 1e5]);
% xlim([-1 17]);
% title('cba18 2 on 18 host');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% subplot(5,4,3)
% errorbar(time,mean(cba18_3_4_phage'),std(cba18_3_4_phage'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e12]);
% yticks([1e1 1e3 1e5 1e7 1e9]);
% xlim([-1 17]);
% title('cba18 3 on 4 phage');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% 
% subplot(5,4,4)
% errorbar(time,mean(cba18_3_4_host'),std(cba18_3_4_host'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e8]);
% yticks([1e1 1e2 1e3 1e4 1e5]);
% xlim([-1 17]);
% title('cba18 3 on 4 host');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% 
% 
% subplot(5,4,5)
% errorbar(time,mean(cba18_3_18_phage'),std(cba18_3_18_phage'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e12]);
% yticks([1e1 1e3 1e5 1e7 1e9]);
% xlim([-1 17]);
% title('cba18 3 on 18 phage');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% 
% subplot(5,4,6)
% errorbar(time,mean(cba18_3_18_host'),std(cba18_3_18_host'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e8]);
% yticks([1e1 1e2 1e3 1e4 1e5]);
% xlim([-1 17]);
% title('cba18 3 on 18 host');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% 
% 
% subplot(5,4,7)
% errorbar(time,mean(cba38_1_38_phage'),std(cba38_1_38_phage'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e12]);
% yticks([1e1 1e3 1e5 1e7 1e9]);
% xlim([-1 17]);
% title('cba38 1 on 38 phage');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% 
% subplot(5,4,8)
% errorbar(time,mean(cba38_1_38_host'),std(cba38_1_38_host'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e8]);
% yticks([1e1 1e2 1e3 1e4 1e5]);
% xlim([-1 17]);
% title('cba38 1 on 38 host');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% 
% 
% subplot(5,4,9)
% errorbar(time,mean(hp1_h100_phage'),std(hp1_h100_phage'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e12]);
% yticks([1e1 1e3 1e5 1e7 1e9]);
% xlim([-1 17]);
% title('hp1 on h100 phage');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% 
% subplot(5,4,10)
% errorbar(time,mean(hp1_h100_host'),std(hp1_h100_host'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e8]);
% yticks([1e0 1e1 1e2 1e3 1e4 1e5]);
% xlim([-1 17]);
% title('hp1 on h100 host');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% 
% subplot(5,4,11)
% errorbar(time,mean(hp1_1315_phage'),std(hp1_1315_phage'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e12]);
% yticks([1e1 1e3 1e5 1e7 1e9]);
% xlim([-1 17]);
% title('hp1 on 13:15 phage');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% 
% subplot(5,4,12)
% errorbar(time,mean(hp1_1315_host'),std(hp1_1315_host'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e8]);
% yticks([1e0 1e1 1e2 1e3 1e4 1e5]);
% xlim([-1 17]);
% title('hp1 on 13:15 host');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% 
% 
% 
% 
% subplot(5,4,13)
% errorbar(time,mean(hs6_h100_phage'),std(hs6_h100_phage'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e12]);
% yticks([1e1 1e3 1e5 1e7 1e9]);
% xlim([-1 17]);
% title('hs6 on h100 phage');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% 
% subplot(5,4,14)
% errorbar(time,mean(hs6_h100_host'),std(hs6_h100_host'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e8]);
% yticks([1e0 1e1 1e2 1e3 1e4 1e5]);
% xlim([-1 17]);
% title('hs6 on h100 host');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% 
% subplot(5,4,15)
% errorbar(time,mean(hs6_1315_phage'),std(hs6_1315_phage'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e12]);
% yticks([1e1 1e3 1e5 1e7 1e9]);
% xlim([-1 17]);
% title('hs6 on 13:15 phage');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% 
% subplot(5,4,16)
% errorbar(time,mean(hs6_1315_host'),std(hs6_1315_host'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e8]);
% yticks([1e0 1e1 1e2 1e3 1e4 1e5]);
% xlim([-1 17]);
% title('hs6 on 13:15 host');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% subplot(5,4,17)
% errorbar(time,mean(cba38_1_18_phage'),std(cba38_1_18_phage'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e12]);
% yticks([1e1 1e3 1e5 1e7 1e9]);
% xlim([-1 17]);
% title('cba 38:1 on 18 phage');
% ylabel('genome copies/ml');
% 
% 
% 
% 
% 
% subplot(5,4,18)
% errorbar(time,mean(cba38_1_18_host'),std(cba38_1_18_host'),LineStyle="none",LineWidth=3);
% xlabel('Time (hr)');
% set(gca,'YScale','log');
% ylim([1e4 1e8]);
% yticks([1e0 1e1 1e2 1e3 1e4 1e5]);
% xlim([-1 17]);
% title('cba 38:1 on 18 host');
% ylabel('genome copies/ml');
% 
