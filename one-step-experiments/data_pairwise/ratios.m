clear all;
clc;
training_time =[];

%% control

cba4_nophage =  1e3*[2.11E+03	1.30E+03	3.00E+03;
2.71E+03	2.39E+03	2.21E+03;
5.50E+04	7.81E+04	7.20E+04];

cba18_nophage = 1e3*[4.89E+03	7.92E+03	9.53E+03;
1.56E+04	9.52E+03	1.26E+04;
3.71E+05	2.39E+05	1.93E+05];

cba38_nophage  = 1e3*[7.81E+03	1.01E+04	7.77E+03;
1.42E+04	1.04E+04	1.29E+04;
1.01E+05	2.02E+05	2.40E+05];

h100_nophage = 1e3*[1.21E+04	1.24E+04	1.18E+04;
7.36E+04	5.34E+04	6.42E+04;
6.52E+05	6.30E+05	9.31E+05];

psa1315_nophage = 1e3*[1.43E+04	1.56E+04	1.84E+04;
2.86E+04	2.31E+04	2.57E+04;
2.95E+05	2.02E+05	3.21E+05];

cba4_nophage_ratio = mean(cba4_nophage(3,:))./mean(cba4_nophage(1,:));
%cba4_nophage_ratio_error = std(cba4_nophage(1,:))/sqrt(3) + std(cba4_nophage(2,:))/sqrt(3);
cba4_nophage_ratio_error = std(cba4_nophage(3,:)./cba4_nophage(1,:))/sqrt(3);


cba18_nophage_ratio = mean(cba18_nophage(3,:))./mean(cba18_nophage(1,:));
cba18_nophage_ratio_error = std(cba18_nophage(3,:)./cba18_nophage(1,:))/sqrt(3);

cba38_nophage_ratio = mean(cba38_nophage(3,:))./mean(cba38_nophage(1,:));
cba38_nophage_ratio_error = std(cba38_nophage(3,:)./cba38_nophage(1,:))/sqrt(3);

h100_nophage_ratio = mean(h100_nophage(3,:))./mean(h100_nophage(1,:));
h100_nophage_ratio_error = std(h100_nophage(3,:)./h100_nophage(1,:))/sqrt(3);

psa1315_nophage_ratio = mean(psa1315_nophage(3,:))./mean(psa1315_nophage(1,:));
psa1315_nophage_ratio_error = std(psa1315_nophage(3,:)./psa1315_nophage(1,:))/sqrt(3);


%% 1
cba18_2_18_host =  1e3*[5.46E+03	6.68E+03	3.78E+03;
9.15E+03	9.78E+03	7.62E+03;
4.01E+04	3.70E+04	4.32E+04];

cba18_2_18_host_ratio = mean(cba18_2_18_host(3,:))./mean(cba18_2_18_host(1,:));
cba18_2_18_host_ratio_error = std(cba18_2_18_host(3,:)./cba18_2_18_host(1,:))/sqrt(3);



%% 2
cba18_3_4_host = 1e3*[3.23E+03	4.18E+03	2.56E+03;
2.20E+03	2.83E+03	1.72E+03;
9.39E+03	3.41E+03	5.98E+03];


cba18_3_4_host_ratio = mean(cba18_3_4_host(3,:))./mean(cba18_3_4_host(1,:));
cba18_3_4_host_ratio_error = std(cba18_3_4_host(3,:)./cba18_3_4_host(1,:))/sqrt(3);


%% 3
cba18_3_18_host = 1e3*[7.87E+03	8.00E+03	6.09E+03;
1.04E+04	1.08E+04	6.60E+03;
1.16E+03	1.36E+03	1.58E+03];


cba18_3_18_host_ratio = mean(cba18_3_18_host(3,:))./mean(cba18_3_18_host(1,:));
cba18_3_18_host_ratio_error = std(cba18_3_18_host(3,:)./cba18_3_18_host(1,:))/sqrt(3);


%% 4

cba38_1_18_host = 1e3*[9.01E+03	8.17E+03	8.13E+03;
1.09E+04	1.35E+04	9.19E+03;
1.18E+03	1.42E+03	1.86E+03];

cba38_1_18_host_ratio = mean(cba38_1_18_host(3,:))./mean(cba38_1_18_host(1,:));
cba38_1_18_host_ratio_error = std(cba38_1_18_host(3,:)./cba38_1_18_host(1,:))/sqrt(3);

%% 5 

cba38_1_38_host = 1e3*[9.01E+03	8.17E+03	8.13E+03;
1.09E+04	1.35E+04	9.19E+03;
1.18E+03	1.42E+03	1.86E+03];

cba38_1_38_host_ratio = mean(cba38_1_38_host(3,:))./mean(cba38_1_38_host(1,:));
cba38_1_38_host_ratio_error = std(cba38_1_38_host(3,:)./cba38_1_38_host(1,:))/sqrt(3);



%% 6

hp1_h100_host = 1e3*[1.03E+04	1.42E+04	1.78E+04;
8.99E+04	8.52E+04	1.30E+05;
9.45E+03	6.21E+03	6.66E+03];



hp1_h100_host_ratio = mean(hp1_h100_host(3,:))./mean(hp1_h100_host(1,:));
hp1_h100_host_ratio_error = std(hp1_h100_host(3,:)./hp1_h100_host(1,:))/sqrt(3);



%% 7
hp1_1315_host = 1e3*[1.48E+04	1.71E+04	2.01E+04
3.89E+04	4.62E+04	4.67E+04
7.19E+04	1.22E+05	8.05E+04];


hp1_1315_host_ratio = mean(hp1_1315_host(3,:))./mean(hp1_1315_host(1,:));
hp1_1315_host_ratio_error = std(hp1_1315_host(3,:)./hp1_1315_host(1,:))/sqrt(3);

%% 8

hs6_h100_host = 1e3*[1.33E+04	1.80E+04	1.69E+04;
7.74E+04	8.50E+04	8.66E+04;
2.35E+05	3.57E+05	2.43E+05];

hs6_h100_host_ratio = mean(hs6_h100_host(3,:))./mean(hs6_h100_host(1,:));
hs6_h100_host_ratio_error = std(hs6_h100_host(3,:)./hs6_h100_host(1,:))/sqrt(3);


%% 9
hs6_1315_host = 1e3*[1.65E+04	2.02E+04	1.99E+04;
3.29E+04	3.08E+04	3.33E+04;
1.76E+05	1.48E+05	1.37E+05];


hs6_1315_host_ratio = mean(hs6_1315_host(3,:))./mean(hs6_1315_host(1,:));
hs6_1315_host_ratio_error = std(hs6_1315_host(3,:)./hs6_1315_host(1,:))/sqrt(3);

%%

uninfected_ratio = [cba18_nophage_ratio,cba4_nophage_ratio,cba18_nophage_ratio,cba38_nophage_ratio,h100_nophage_ratio,psa1315_nophage_ratio,h100_nophage_ratio,psa1315_nophage_ratio];
uninfected_ratio_error = [cba18_nophage_ratio_error,cba4_nophage_ratio_error,cba18_nophage_ratio_error,cba38_nophage_ratio_error,h100_nophage_ratio_error,psa1315_nophage_ratio_error,h100_nophage_ratio_error,psa1315_nophage_ratio_error];

infected_ratio = [cba18_2_18_host_ratio, cba18_3_4_host_ratio,cba18_3_18_host_ratio, cba38_1_38_host_ratio, hp1_h100_host_ratio,hp1_1315_host_ratio,hs6_h100_host_ratio,hs6_1315_host_ratio];
infected_ratio_error = [cba18_2_18_host_ratio_error, cba18_3_4_host_ratio_error,cba18_3_18_host_ratio_error, cba38_1_38_host_ratio_error, hp1_h100_host_ratio_error,hp1_1315_host_ratio_error,hs6_h100_host_ratio_error,hs6_1315_host_ratio_error];
%%


% uninfected_ratio = [cba4_nophage_ratio,cba38_nophage_ratio,cba18_nophage_ratio,cba18_nophage_ratio,psa1315_nophage_ratio,psa1315_nophage_ratio,h100_nophage_ratio,h100_nophage_ratio];
% uninfected_ratio_error = [cba4_nophage_ratio_error,cba38_nophage_ratio_error,cba18_nophage_ratio_error,cba18_nophage_ratio_error,psa1315_nophage_ratio_error,psa1315_nophage_ratio_error,h100_nophage_ratio_error,h100_nophage_ratio_error];
% 
% infected_ratio = [cba18_3_4_host_ratio,cba38_1_38_host_ratio,cba18_2_18_host_ratio,cba18_3_18_host_ratio, hp1_1315_host_ratio, hs6_1315_host_ratio, hp1_h100_host_ratio, hs6_h100_host_ratio];
% infected_ratio_error = [cba18_3_4_host_ratio_error,cba38_1_38_host_ratio_error,cba18_2_18_host_ratio_error,cba18_3_18_host_ratio_error,hp1_1315_host_ratio_error, hs6_1315_host_ratio_error, hp1_h100_host_ratio_error, hs6_h100_host_ratio_error];
% 
% 
% 
% %infected_ratio = [cba18_3_4_host_ratio,cba38_1_38_host_ratio,cba18_2_18_host_ratio,cba18_3_18_host_ratio,hs6_h100_host_ratio,hs6_1315_host_ratio,hp1_h100_host_ratio,hp1_1315_host_ratio];
% %infected_ratio_error = [cba18_3_4_host_ratio_error,cba38_1_38_host_ratio_error,cba18_2_18_host_ratio_error,cba18_3_18_host_ratio_error,hs6_h100_host_ratio_error,hs6_1315_host_ratio_error,hp1_h100_host_ratio_error,hp1_1315_host_ratio_error];

%% figure
x_axis = 1:2:16;
gap = 0.1
figure(1)
errorbar(x_axis-gap,uninfected_ratio,uninfected_ratio_error,"Marker","o", ...
    "MarkerEdgeColor","k","MarkerFaceColor","none","MarkerSize",15, ...
    "LineStyle","none","Color","k");

hold on;

errorbar(x_axis,infected_ratio,infected_ratio_error,"Marker","o", ...
    "MarkerEdgeColor","k","MarkerFaceColor","k","MarkerSize",15, ...
    "LineStyle","none","Color","k");

set(gca, 'FontSize', 20);
set(gca, 'YScale', 'log');
set(gca,'fontname','times');
xlim([0 16]);
xticks(x_axis)



%%



cba18_2_18_host = 1e3*[5.46E+03	6.68E+03	3.78E+03;
9.15E+03	9.78E+03	7.62E+03;
4.01E+04	3.70E+04	4.32E+04];

cba18_2_18_phage = 1e3*[1.73E+03	1.89E+03	2.53E+03;
1.56E+04	1.10E+04	3.18E+04;
1.16E+08	1.26E+08	1.14E+08];



cba18_3_4_host = 1e3*[3.23E+03	4.18E+03	2.56E+03;
2.20E+03	2.83E+03	1.72E+03;
9.39E+03	3.41E+03	5.98E+03];



cba18_3_4_phage = 1e3*[4.46E+02	4.37E+02	4.40E+02;
5.33E+02	5.30E+02	5.15E+02;
3.69E+05	3.89E+05	4.54E+05];


cba18_3_18_host = 1e3*[7.87E+03	8.00E+03	6.09E+03;
1.04E+04	1.08E+04	6.60E+03;
1.16E+03	1.36E+03	1.58E+03];


cba18_3_18_phage = 1e3*[4.54E+02	4.71E+02	5.45E+02;
3.78E+03	2.70E+03	2.99E+03;
1.57E+06	1.66E+06	1.37E+06];


cba38_1_38_host = 1e3*[9.01E+03	8.17E+03	8.13E+03;
1.09E+04	1.35E+04	9.19E+03
1.18E+03	1.42E+03	1.86E+03];


cba38_1_38_phage = 1e3*[2.85E+02	2.73E+02	2.90E+02;
2.65E+03	2.48E+03	2.23E+03;
2.81E+05	3.08E+05	3.20E+05];


hp1_h100_host = 1e3*[1.03E+04	1.42E+04	1.78E+04;
8.99E+04	8.52E+04	1.30E+05;
9.45E+03	6.21E+03	6.66E+03];

hp1_h100_phage = 1e3*[3.03E+02	1.94E+02	4.94E+02;
1.38E+04	7.82E+03	1.07E+04;
8.77E+07	7.95E+07	8.48E+07];


hp1_1315_host = 1e3*[1.48E+04	1.71E+04	2.01E+04
3.89E+04	4.62E+04	4.67E+04
7.19E+04	1.22E+05	8.05E+04];


hp1_1315_phage = 1e3*[4.44E+02	2.82E+02	4.34E+02;
1.07E+04	9.97E+03	1.61E+04;
5.53E+07	5.50E+07	4.84E+07];


hs6_h100_host =1e3*[1.33E+04	1.80E+04	1.69E+04;
7.74E+04	8.50E+04	8.66E+04;
2.35E+05	3.57E+05	2.43E+05];

hs6_h100_phage = 1e3*[7.79E+02	5.53E+02	7.96E+02;
3.95E+04	2.52E+04	1.86E+04;
1.12E+07	1.10E+07	9.59E+06];


hs6_1315_host =  1e3*[1.65E+04	2.02E+04	1.99E+04;
3.29E+04	3.08E+04	3.33E+04;
1.76E+05	1.48E+05	1.37E+05];

hs6_1315_phage = 1e3*[7.06E+02	7.25E+02	6.49E+02;
2.58E+04	2.93E+04	2.65E+04;
1.07E+07	1.25E+07	1.32E+07];

%% 1

qprc_host = cba18_2_18_host ;
qprc_virus = cba18_2_18_phage;



load('./../data_2024/CBA18-2_18_2024.mat');
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
virus_den = y_series_opt(NE+3,:);
clear y0

cba18_2_18_host_ratio_model = host_den(end)./host_den(1);
theta

training_time(end+1) =  time_free_phages(end);
%% 2
qprc_host = cba18_3_4_host ;
qprc_virus = cba18_3_4_phage;
load('./../data_2024/CBA18-3_4_2024.mat');
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
virus_den = y_series_opt(NE+3,:);
clear y0

cba18_3_4_host_ratio_model = host_den(end)./host_den(1);
theta

training_time(end+1) =  time_free_phages(end)/60;
%% 3

qprc_host = cba18_3_18_host ;
qprc_virus = cba18_3_18_phage;

load('./../data_2024/CBA18-3_18_2024.mat');


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
virus_den = y_series_opt(NE+3,:);
clear y0

cba18_3_18_host_ratio_model = host_den(end)./host_den(1);
theta
training_time(end+1) =  time_free_phages(end)*60;

%%  4

qprc_host = cba38_1_38_host ;
qprc_virus = cba38_1_38_phage;

load('./../data_2024/CBA38-1_38.mat');

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
virus_den = y_series_opt(NE+3,:);
clear y0

cba38_1_38_host_ratio_model = host_den(end)./host_den(1);
theta
training_time(end+1) =  time_free_phages(end)*60;


%%  5

qprc_host = hp1_h100_host;
qprc_virus = hp1_h100_phage;

load('./../data_2024/HP1_H100_2024.mat');

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
virus_den = y_series_opt(NE+3,:);
clear y0



hp1_h100_host_ratio_model = host_den(end)./host_den(1);
theta
training_time(end+1) =  time_free_phages(end);

%% 6

qprc_host = hp1_1315_host;
qprc_virus = hp1_1315_phage;

load('./../data_2024/HP1_13-15_2024.mat');

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
virus_den = y_series_opt(NE+3,:);
clear y0

hp1_1315_host_ratio_model = host_den(end)./host_den(1);
theta

training_time(end+1) =  time_free_phages(end)*60;


%%  7

qprc_host = hs6_h100_host;
qprc_virus = hs6_h100_phage;

load('./../data_2024/HS6_H100_2024.mat');

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
virus_den = y_series_opt(NE+3,:);
clear y0

hp6_h100_host_ratio_model = host_den(end)./host_den(1);
theta

training_time(end+1) =  time_free_phages(end);

%%  8


qprc_host = hs6_1315_host;
qprc_virus = hs6_1315_phage;

load('./../data_2024/HS6_13-15_2024.mat');
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
virus_den = y_series_opt(NE+3,:);
clear y0

hp6_1315_host_ratio_model = host_den(end)/host_den(1);

theta
training_time(end+1) =  time_free_phages(end);

%%

%model_ratio = [cba18_3_4_host_ratio_model,cba38_1_38_host_ratio_model,cba18_2_18_host_ratio_model,cba18_3_18_host_ratio_model,hp6_h100_host_ratio_model, hp6_1315_host_ratio_model, hp1_h100_host_ratio_model , hp1_1315_host_ratio_model];
model_ratio = [cba18_2_18_host_ratio_model,cba18_3_4_host_ratio_model, cba18_3_18_host_ratio_model,cba38_1_38_host_ratio_model,hp1_h100_host_ratio_model,hp1_1315_host_ratio_model,hp6_h100_host_ratio_modelfig = gcf;,hp6_1315_host_ratio_model];
limit_of_detection = 1e-3;

model_ratio(model_ratio<limit_of_detection) = 1e-3;
%%
figure(1)
plot(x_axis+gap,model_ratio,"Marker","diamond", ...
    "MarkerEdgeColor","k","MarkerFaceColor",[0.8,0.8,0.8],"MarkerSize",15, ...
    "LineStyle","none","Color","k");
%yline(limit_of_detection,"LineWidth",3,"LineStyle","--")
yticks([1e-3,1e-2,1e-1,1e0,1e1,1e2]);
yticklabels({'<10^{-3}','10^{-2}','10^{-1}','1','10','10^2'})

ylim([1e-3,1e2]);
yline(1,"LineWidth",3,"LineStyle","-.")

%xticklabels({'\phi18:3 – CBA 4  ','\phi38:1 – CBA 38  ','\phi18:2 – CBA 18  ','\phi18:3 – CBA 18  ','PSA-HP1 – PSA H100  ','PSA-HP1 – PSA 13-15  ','PSA-HS6 – PSA H100  ','PSA-HS6 – PSA 13-15'})
xticklabels({'\phi18:2 – CBA 18  ','\phi18:3 – CBA 4  ','\phi18:3 – CBA 18  ','\phi38:1 – CBA 38  ','PSA-HP1 – PSA H100  ','PSA-HP1 – PSA 13-15  ','PSA-HS6 – PSA H100  ','PSA-HS6 – PSA 13-15'})

xtickangle(90)
ylabel("Host_{final}/Host_{initial}");
legend('Uninfected (data)','Infected (data)','Pairwise SEIV model')
%set(gca,'TickLength',[0.1 0.1])

%line([0 0.5], [1e2 1e2], '-k', 'LineWidth',1);

%% figure 2

