clc;
clear all;

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



figure(2)
options = odeset('NonNegative',1);

%%



qprc_host = cba18_3_4_host ;
qprc_virus = cba18_3_4_phage;
load('./data_2024/CBA18-3_4_2024.mat');
name = string(labels.phage)+'_'+string(labels.host);
load('parameters.mat','pars');

data.xdata{1} = time_free_phages;
data.ydata{1} = free_phages;



NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(qprc_host(1,:));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qprc_virus(1,:));

tvec = 0:0.1:15.75;

NE = round(theta_optimized(5));
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0,options);

y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));
clear y0

subplot(2,4,1)
plot(time_opt,host_den,'LineWidth',3,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('\phi18:3—CBA 4');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')




%%

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
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0,options);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));

subplot(2,4,2)
plot(time_opt,host_den,'LineWidth',3,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('\phi38:1—CBA 38');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')


%%

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
options = odeset('NonNegative',1);
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0,options);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));
subplot(2,4,3)
plot(time_opt,host_den,'LineWidth',3,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('\phi18:2—CBA 18');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')


%%
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
tvec_short = 0:0.1:15.75;

theta
NE = round(theta_optimized(5));
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0,options);


y_series_opt = deval(solved,tvec);

%y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));

virus_den = y_series_opt(NE+3,:);

%subplot(2,4,4)
%plot(time_opt(1:10:end),host_den(1:10:end),'.','LineWidth',3,'Color',[0.8,0.8,0.8],'LineStyle','none');
plot(tvec,host_den,'.','LineWidth',3,'Color',[0.8,0.8,0.8],'LineStyle','none');

set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('\phi18:3—CBA 18');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')


%%

subplot(2,2,1)
plot(tvec,host_den,'.');
xlabel('Time (hrs)');
ylabel('Host density (genome copies/ml)');
xlabel('Time (hrs)');
ylabel('Host density (genome copies/ml)');
title('Host linear scale');
subplot(2,2,2)
plot(tvec,host_den,'.');
set(gca, 'YScale', 'log');
xlabel('Time (hrs)');
ylabel('Host density (genome copies/ml)');
xlabel('Time (hrs)');
ylabel('Host density (genome copies/ml)');
title('Host log scale');
subplot(2,2,3)
plot(tvec,virus_den,'.');
xlabel('Time (hrs)');
ylabel('Phage density (genome copies/ml)');
xlabel('Time (hrs)');
ylabel('Host density (genome copies/ml)');
title('Phage linear scale');
subplot(2,2,4)
plot(tvec,virus_den,'.');
set(gca, 'YScale', 'log');
xlabel('Time (hrs)');
ylabel('Phage density (genome copies/ml)');
xlabel('Time (hrs)');
ylabel('Host density (genome copies/ml)');
title('Phage log scale');


%%

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
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0,options);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));

subplot(2,4,5)
plot(time_opt,host_den,'LineWidth',3,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('PSA HS6—PSA H100');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')


%%


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
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0,options);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));
subplot(2,4,6)
plot(time_opt,host_den,'LineWidth',3,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('PSA HS6—PSA 13:15');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')


%%

qprc_host = hp1_h100_host;
qprc_virus = hp1_h100_phage;

load('./data_2024/HP1_H100_2024.mat');

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');

theta = theta_optimized;


NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(qprc_host(1,:));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qprc_virus(1,:));

tvec = 0:0.1:15.75;

NE = round(theta_optimized(5));
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0,options);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));
subplot(2,4,7)
plot(time_opt,host_den,'LineWidth',3,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('PSA HP1—PSA H100');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')



%%

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
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0,options);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));

subplot(2,4,8)
plot(time_opt,host_den,'.','LineWidth',3,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('PSA HP1—PSA 13:15');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')