clc;
clear all;
addpath(genpath(pwd)); % add current directory to path (mcmcstat is included)
addpath('./simulator/');
addpath('./mcmcstat/');

%%

cba18_2_18_host = 1e3*[4.30E+02	2.66E+02	4.71E+02;
2.90E+03	2.27E+03	2.71E+03];

cba18_2_18_phage = 1e3*[1.50E+04	1.50E+04	1.37E+04;
6.01E+08	5.57E+08	5.58E+08];

cba18_3_4_host = 1e3*[1.04E+02	1.14E+02	2.28E+02;
5.37E+02	5.26E+02	5.16E+02];

cba18_3_4_phage = 1e3*[1.11E+04	1.30E+04	1.29E+04;
1.97E+07	1.91E+07	1.90E+07];

cba18_3_18_host = 1e3*[3.32E+02	2.26E+02	3.80E+02;
1.45E+02	1.02E+02	7.51E+01];

cba18_3_18_phage = 1e3*[1.08E+04	1.12E+04	1.22E+04;
6.36E+07	5.93E+07	4.98E+07];

cba38_1_38_host = 1e3*[1.07E+02	1.15E+02	1.00E+02;
1.28E+02	8.30E+01	2.68E+02];

cba38_1_38_phage = 1e3*[3.86E+05	4.23E+05	4.22E+05;
9.79E+06	9.57E+06	9.60E+06];


cba38_1_18_host = 1e3*[2.35E+02	3.26E+02	3.33E+02;
3.87E+03	3.29E+03	3.80E+03];

cba38_1_18_phage = 1e3*[4.32E+05	4.75E+05	4.76E+05; %this had a mistake
4.48E+05	4.46E+05	4.00E+05];

hp1_h100_host = 1e3*[8.79E+01	5.59E+01	2.70E+01;
3.05E+01	3.60E+00	2.01E+01];

hp1_h100_phage = 1e3*[2.55E+03	3.41E+03	3.51E+03;
8.25E+07	7.98E+07	7.45E+07];

hp1_1315_host = 1e3*[1.58E+02	9.28E+01	2.59E+02;
8.97E+00	7.95E+00	1.01E+01];

hp1_1315_phage = 1e3*[4.22E+03	3.60E+03	3.37E+03;
5.54E+07	5.31E+07	5.46E+07];

hs6_h100_host = 1e3*[1.16E+01	6.26E+01	6.66E+01;
9.89E+02	2.36E+02	9.40E+01];

hs6_h100_phage = 1e3*[1.05E+04	6.57E+06	1.09E+04;
2.08E+08	1.68E+08	1.85E+08];

hs6_1315_host = 1e3*[9.02E+01	4.22E+02	1.03E+02;
6.25E+02	7.06E+02	1.30E+03];

hs6_1315_phage = 1e3*[1.02E+04	1.07E+04	1.21E+04;
1.03E+08	9.02E+07	9.36E+07];

%% 1

qprc_host = cba18_2_18_host ;
qprc_virus = cba18_2_18_phage;



load('./data_2024/CBA18-2_18_2024.mat');
name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');


NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(qprc_host(1,:));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qprc_virus(1,:));

tvec = 0:0.1:15.75;

NE = round(theta_optimized(5));
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));
virus_den = y_series_opt(NE+3,:);
clear y0




subplot(4,4,1)
plot(time_opt,virus_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virions/ml')
title('CBA18-2 on CBA 18');
set(gca,'FontSize',20);
xlim([0 15.75]);
yticks([1e3 1e5 1e7 1e9 1e11]);
errorbar([0 945]/60,mean(cba18_2_18_phage'),std(cba18_2_18_phage'),  'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','b','LineStyle','none', 'LineWidth',3, 'Color','b');
ylim([1e6 1e12]);



subplot(4,4,2)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('hosts/ml');
set(gca,'FontSize',20);
title('CBA18-2 on CBA 18');
xlim([0 15.75]);
ylim([1e-1 1e9]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);
errorbar([0 945]/60,mean(cba18_2_18_host'),std(cba18_2_18_host'),  'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','LineStyle','none', 'LineWidth',3,'Color','r');



%% 2

qprc_host = cba18_3_4_host ;
qprc_virus = cba18_3_4_phage;

load('./data_2024/CBA18-3_4_2024.mat');

name = string(labels.phage)+'_'+string(labels.host);
load('parameters.mat','pars');


format short g



data.xdata{1} = time_free_phages;
data.ydata{1} = free_phages;


NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(qprc_host(1,:));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qprc_virus(1,:));

tvec = 0:0.1:15.75;

NE = round(theta_optimized(5));
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));
virus_den = y_series_opt(NE+3,:);
clear y0




subplot(4,4,3)
plot(time_opt,virus_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virions/ml')
title('CBA18-3 on CBA 4');
set(gca,'FontSize',20);
xlim([0 15.75]);
yticks([1e3 1e5 1e7 1e9 1e11]);
errorbar([0 945]/60,mean(qprc_virus'),std(qprc_virus'),  'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','b','LineStyle','none', 'LineWidth',3, 'Color','b');
ylim([1e6 1e12]);



subplot(4,4,4)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('hosts/ml');
set(gca,'FontSize',20);
title('CBA18-3 on CBA 4');
xlim([0 15.75]);
ylim([1e-1 1e9]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);
errorbar([0 945]/60,mean(qprc_host'),std(qprc_host'),  'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','LineStyle','none', 'LineWidth',3,'Color','r');


%% 3

qprc_host = cba18_3_18_host ;
qprc_virus = cba18_3_18_phage;

load('./data_2024/CBA18-3_18_2024.mat');


name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');

NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(qprc_host(1,:));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qprc_virus(1,:));

tvec = 0:0.1:15.75;

NE = round(theta_optimized(5));
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));
virus_den = y_series_opt(NE+3,:);
clear y0




subplot(4,4,5)
plot(time_opt,virus_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virions/ml')
title('CBA18-3 on CBA 18');
set(gca,'FontSize',20);
xlim([0 15.75]);
yticks([1e3 1e5 1e7 1e9 1e11]);
errorbar([0 945]/60,mean(qprc_virus'),std(qprc_virus'),  'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','b','LineStyle','none', 'LineWidth',3, 'Color','b');
ylim([1e6 1e12]);



subplot(4,4,6)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('hosts/ml');
set(gca,'FontSize',20);
title('CBA18-3 on CBA 18');
xlim([0 15.75]);
ylim([1e-1 1e9]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);
errorbar([0 945]/60,mean(qprc_host'),std(qprc_host'),  'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','LineStyle','none', 'LineWidth',3,'Color','r');




%%  4

qprc_host = cba38_1_38_host ;
qprc_virus = cba38_1_38_phage;

load('./data_2024/CBA38-1_38.mat');

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');

NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(qprc_host(1,:));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qprc_virus(1,:));

tvec = 0:0.1:15.75;

NE = round(theta_optimized(5));
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));
virus_den = y_series_opt(NE+3,:);
clear y0




subplot(4,4,7)
plot(time_opt,virus_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virions/ml')
title('CBA38-1 on CBA 38');
set(gca,'FontSize',20);
xlim([0 15.75]);
yticks([1e3 1e5 1e7 1e9 1e11]);
errorbar([0 945]/60,mean(qprc_virus'),std(qprc_virus'),  'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','b','LineStyle','none', 'LineWidth',3, 'Color','b');
ylim([1e6 1e12]);



subplot(4,4,8)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('hosts/ml');
set(gca,'FontSize',20);
title('CBA38-1 on CBA 38');
xlim([0 15.75]);
ylim([1e-1 1e9]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);
errorbar([0 945]/60,mean(qprc_host'),std(qprc_host'),  'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','LineStyle','none', 'LineWidth',3,'Color','r');



%%  5

qprc_host = hp1_h100_host;
qprc_virus = hp1_h100_phage;

load('./data_2024/HP1_H100_2024.mat');

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');




NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(qprc_host(1,:));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qprc_virus(1,:));

tvec = 0:0.1:15.75;

NE = round(theta_optimized(5));
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));
virus_den = y_series_opt(NE+3,:);
clear y0




subplot(4,4,9)
plot(time_opt,virus_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virions/ml')
title('HP1 on H100');
set(gca,'FontSize',20);
xlim([0 15.75]);
yticks([1e3 1e5 1e7 1e9 1e11]);
errorbar([0 945]/60,mean(qprc_virus'),std(qprc_virus'),  'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','b','LineStyle','none', 'LineWidth',3, 'Color','b');
ylim([1e6 1e12]);



subplot(4,4,10)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('hosts/ml');
set(gca,'FontSize',20);
title('HP1 on H100');
xlim([0 15.75]);
ylim([1e-1 1e9]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);
errorbar([0 945]/60,mean(qprc_host'),std(qprc_host'),  'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','LineStyle','none', 'LineWidth',3,'Color','r');




%% 6

qprc_host = hp1_1315_host;
qprc_virus = hp1_1315_phage;

load('./data_2024/HP1_13-15_2024.mat');

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');

y0(1) = mean(qprc_host(1,:));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qprc_virus(1,:));

tvec = 0:0.1:15.75;

NE = round(theta_optimized(5));
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));
virus_den = y_series_opt(NE+3,:);
clear y0




subplot(4,4,11)
plot(time_opt,virus_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virions/ml')
title('HP1 on 1315');
set(gca,'FontSize',20);
xlim([0 15.75]);
yticks([1e3 1e5 1e7 1e9 1e11]);
errorbar([0 945]/60,mean(qprc_virus'),std(qprc_virus'),  'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','b','LineStyle','none', 'LineWidth',3, 'Color','b');
ylim([1e6 1e12]);



subplot(4,4,12)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('hosts/ml');
set(gca,'FontSize',20);
title('HP1 on 1315');
xlim([0 15.75]);
ylim([1e-1 1e9]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);
errorbar([0 945]/60,mean(qprc_host'),std(qprc_host'),  'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','LineStyle','none', 'LineWidth',3,'Color','r');


%%  7

qprc_host = hs6_h100_host;
qprc_virus = hs6_h100_phage;

load('./data_2024/HS6_H100_2024.mat');

y0(1) = mean(qprc_host(1,:));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qprc_virus(1,:));
tvec = 0:0.1:15.75;

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');



NE = round(theta_optimized(5));
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));
virus_den = y_series_opt(NE+3,:);
clear y0




subplot(4,4,13)
plot(time_opt,virus_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virions/ml')
title('HS6 on H100');
set(gca,'FontSize',20);
xlim([0 15.75]);
yticks([1e3 1e5 1e7 1e9 1e11]);
errorbar([0 945]/60,mean(qprc_virus'),std(qprc_virus'),  'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','b','LineStyle','none', 'LineWidth',3, 'Color','b');
ylim([1e6 1e12]);


subplot(4,4,14)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('hosts/ml');
set(gca,'FontSize',20);
title('HS6 on HS100');
xlim([0 15.75]);
ylim([1e-1 1e9]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);
errorbar([0 945]/60,mean(qprc_host'),std(qprc_host'),  'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','LineStyle','none', 'LineWidth',3,'Color','r');




%%  8


qprc_host = hs6_1315_host;
qprc_virus = hs6_1315_phage;

load('./data_2024/HS6_13-15_2024.mat');
y0(1) = mean(qprc_host(1,:));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qprc_virus(1,:));

tvec = 0:0.1:15.75;

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');



NE = round(theta_optimized(5));
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));
virus_den = y_series_opt(NE+3,:);
clear y0




subplot(4,4,15)
plot(time_opt,virus_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virions/ml')
title('HS6 on 1315');
set(gca,'FontSize',20);
xlim([0 15.75]);
yticks([1e3 1e5 1e7 1e9 1e11]);
errorbar([0 945]/60,mean(qprc_virus'),std(qprc_virus'),  'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','b','LineStyle','none', 'LineWidth',3, 'Color','b');
ylim([1e6 1e12]);

subplot(4,4,16)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('hosts/ml');
set(gca,'FontSize',20);
title('HS6 on 1315');
xlim([0 15.75]);
ylim([1e-1 1e9]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);
errorbar([0 945]/60,mean(qprc_host'),std(qprc_host'),  'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','LineStyle','none', 'LineWidth',3,'Color','r');



% Display the legend for this subplot
%legend('Virus Density Baseline','Ineffective viral adsorption','Ineffective viral infection','Virus Data','AutoUpdate','off');
%legend('show');
