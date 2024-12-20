clear all;
clc;
training_time =[];

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

cba4_nophage_ratio = mean(cba4_nophage(2,:))./mean(cba4_nophage(1,:));
%cba4_nophage_ratio_error = std(cba4_nophage(1,:))/sqrt(3) + std(cba4_nophage(2,:))/sqrt(3);
cba4_nophage_ratio_error = std(cba4_nophage(2,:)./cba4_nophage(1,:))/sqrt(3);


cba18_nophage_ratio = mean(cba18_nophage(2,:))./mean(cba18_nophage(1,:));
cba18_nophage_ratio_error = std(cba18_nophage(2,:)./cba18_nophage(1,:))/sqrt(3);

cba38_nophage_ratio = mean(cba38_nophage(2,:))./mean(cba38_nophage(1,:));
cba38_nophage_ratio_error = std(cba38_nophage(2,:)./cba38_nophage(1,:))/sqrt(3);

h100_nophage_ratio = mean(h100_nophage(2,:))./mean(h100_nophage(1,:));
h100_nophage_ratio_error = std(h100_nophage(2,:)./h100_nophage(1,:))/sqrt(3);

psa1315_nophage_ratio = mean(psa1315_nophage(2,:))./mean(psa1315_nophage(1,:));
psa1315_nophage_ratio_error = std(psa1315_nophage(2,:)./psa1315_nophage(1,:))/sqrt(3);


%% 1
cba18_2_18_host = 1e3*[4.30E+02	2.66E+02	4.71E+02;
2.90E+03	2.27E+03	2.71E+03];

cba18_2_18_host_ratio = mean(cba18_2_18_host(2,:))./mean(cba18_2_18_host(1,:));
cba18_2_18_host_ratio_error = std(cba18_2_18_host(2,:)./cba18_2_18_host(1,:))/sqrt(3);



%% 2
cba18_3_4_host = 1e3*[1.04E+02	1.14E+02	2.28E+02;
5.37E+02	5.26E+02	5.16E+02];

cba18_3_4_host_ratio = mean(cba18_3_4_host(2,:))./mean(cba18_3_4_host(1,:));
cba18_3_4_host_ratio_error = std(cba18_3_4_host(2,:)./cba18_3_4_host(1,:))/sqrt(3);


%% 3
cba18_3_18_host = 1e3*[3.32E+02	2.26E+02	3.80E+02;
1.45E+02	1.02E+02	7.51E+01];

cba18_3_18_host_ratio = mean(cba18_3_18_host(2,:))./mean(cba18_3_18_host(1,:));
cba18_3_18_host_ratio_error = std(cba18_3_18_host(2,:)./cba18_3_18_host(1,:))/sqrt(3);


%% 4

cba38_1_18_host = 1e3*[2.35E+02	3.26E+02	3.33E+02;
3.87E+03	3.29E+03	3.80E+03];

cba38_1_18_host_ratio = mean(cba38_1_18_host(2,:))./mean(cba38_1_18_host(1,:));
cba38_1_18_host_ratio_error = std(cba38_1_18_host(2,:)./cba38_1_18_host(1,:))/sqrt(3);

%% 5 

cba38_1_38_host = 1e3*[1.07E+02	1.15E+02	1.00E+02;
1.28E+02	8.30E+01	2.68E+02];

cba38_1_38_host_ratio = mean(cba38_1_38_host(2,:))./mean(cba38_1_38_host(1,:));
cba38_1_38_host_ratio_error = std(cba38_1_38_host(2,:)./cba38_1_38_host(1,:))/sqrt(3);


%% 6

hp1_h100_host = 1e3*[8.79E+01	5.59E+01	2.70E+01;
3.05E+01	3.60E+00	2.01E+01];

hp1_h100_host_ratio = mean(hp1_h100_host(2,:))./mean(hp1_h100_host(1,:));
hp1_h100_host_ratio_error = std(hp1_h100_host(2,:)./hp1_h100_host(1,:))/sqrt(3);



%% 7
hp1_1315_host = 1e3*[1.58E+02	9.28E+01	2.59E+02;
8.97E+00	7.95E+00	1.01E+01];

hp1_1315_host_ratio = mean(hp1_1315_host(2,:))./mean(hp1_1315_host(1,:));
hp1_1315_host_ratio_error = std(hp1_1315_host(2,:)./hp1_1315_host(1,:))/sqrt(3);

%% 8

hs6_h100_host = 1e3*[1.16E+01	6.26E+01	6.66E+01;
9.89E+02	2.36E+02	9.40E+01];

hs6_h100_host_ratio = mean(hs6_h100_host(2,:))./mean(hs6_h100_host(1,:));
hs6_h100_host_ratio_error = std(hs6_h100_host(2,:)./hs6_h100_host(1,:))/sqrt(3);


%% 9
hs6_1315_host = 1e3*[9.02E+01	4.22E+02	1.03E+02;
6.25E+02	7.06E+02	1.30E+03];

hs6_1315_host_ratio = mean(hs6_1315_host(2,:))./mean(hs6_1315_host(1,:));
hs6_1315_host_ratio_error = std(hs6_1315_host(2,:)./hs6_1315_host(1,:))/sqrt(3);


uninfected_ratio = [cba4_nophage_ratio,cba38_nophage_ratio,cba18_nophage_ratio,cba18_nophage_ratio,psa1315_nophage_ratio,psa1315_nophage_ratio,h100_nophage_ratio,h100_nophage_ratio];
uninfected_ratio_error = [cba4_nophage_ratio_error,cba38_nophage_ratio_error,cba18_nophage_ratio_error,cba18_nophage_ratio_error,psa1315_nophage_ratio_error,psa1315_nophage_ratio_error,h100_nophage_ratio_error,h100_nophage_ratio_error];

infected_ratio = [cba18_3_4_host_ratio,cba38_1_38_host_ratio,cba18_2_18_host_ratio,cba18_3_18_host_ratio,hs6_h100_host_ratio,hs6_1315_host_ratio,hp1_h100_host_ratio,hp1_1315_host_ratio];
infected_ratio_error = [cba18_3_4_host_ratio_error,cba38_1_38_host_ratio_error,cba18_2_18_host_ratio_error,cba18_3_18_host_ratio_error,hs6_h100_host_ratio_error,hs6_1315_host_ratio_error,hp1_h100_host_ratio_error,hp1_1315_host_ratio_error];


%% figure
x_axis = 1:2:16;
gap = 0.1
figure(1)
errorbar(x_axis-gap,uninfected_ratio,uninfected_ratio_error,"Marker","o", ...
    "MarkerEdgeColor","k","MarkerFaceColor","none","MarkerSize",10, ...
    "LineStyle","none","Color","k");

hold on;

errorbar(x_axis,infected_ratio,infected_ratio_error,"Marker","o", ...
    "MarkerEdgeColor","k","MarkerFaceColor","k","MarkerSize",10, ...
    "LineStyle","none","Color","k");

set(gca, 'FontSize', 20);
set(gca, 'YScale', 'log');
set(gca,'fontname','times');
xlim([0 16]);
xticks(x_axis)



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
virus_den = y_series_opt(NE+3,:);
clear y0

cba18_3_4_host_ratio_model = host_den(end)./host_den(1);
theta

training_time(end+1) =  time_free_phages(end)/60;
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
virus_den = y_series_opt(NE+3,:);
clear y0

cba38_1_38_host_ratio_model = host_den(end)./host_den(1);
theta
training_time(end+1) =  time_free_phages(end)*60;


%%  5

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
virus_den = y_series_opt(NE+3,:);
clear y0



hp1_h100_host_ratio_model = host_den(end)./host_den(1);
theta
training_time(end+1) =  time_free_phages(end);

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
virus_den = y_series_opt(NE+3,:);
clear y0

hp6_h100_host_ratio_model = host_den(end)./host_den(1);
theta

training_time(end+1) =  time_free_phages(end);

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

model_ratio = [cba18_3_4_host_ratio_model,cba38_1_38_host_ratio_model,cba18_2_18_host_ratio_model,cba18_3_18_host_ratio_model,hp6_h100_host_ratio_model, hp6_1315_host_ratio_model, hp1_h100_host_ratio_model , hp1_1315_host_ratio_model];

limit_of_detection = 1e-7;

model_ratio(model_ratio<limit_of_detection) = 1e-7;
%%
figure(1)
plot(x_axis+gap,model_ratio,"Marker","diamond", ...
    "MarkerEdgeColor","k","MarkerFaceColor",[0.8,0.8,0.8],"MarkerSize",10, ...
    "LineStyle","none","Color","k");
%yline(limit_of_detection,"LineWidth",3,"LineStyle","--")
yticks([1e-7 1e-6 1e-5 1e-4 1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4]);
yticklabels({'<10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','1','10','10^2','10^3','10^4'})

ylim([1e-7,1e4]);
yline(1,"LineWidth",3,"LineStyle","-.")

xticklabels({'\phi18:3—CBA 4  ','\phi38:1—CBA 38  ','\phi18:2—CBA 18  ','\phi18:3—CBA 18  ','PSA HS6—PSA H100  ','PSA HS6—PSA 13-15  ','PSA HP1—PSA H100  ','PSA HP1—PSA 13-15  '})
xtickangle(90)
ylabel("Host_{final}/Host_{initial}");
legend('Uninfected (data)','Infected (data)','Pairwise model')
%set(gca,'TickLength',[0.1 0.1])

%line([0 0.5], [1e2 1e2], '-k', 'LineWidth',1);

%% figure 2



figure(2)





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
xline(training_time(1),'k--'); hold on;
plot(time_opt,host_den,'LineWidth',5,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('\phi18:3--CBA 4','Interpreter','latex');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')
box on;




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
xline(training_time(2),'k--'); hold on;
plot(time_opt,host_den,'LineWidth',5,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('\phi38:1--CBA 38','Interpreter','latex');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')
box on;



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
xline(training_time(3),'k--'); hold on;
plot(time_opt,host_den,'LineWidth',5,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('\phi18:2--CBA 18','Interpreter','latex');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')
box on;


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
solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tvec',y0,options);
y_series_opt = solved.y;
time_opt=solved.x;

host_den = sum(y_series_opt(1:NE+2,:));
subplot(2,4,4)
xline(training_time(4),'k--'); hold on;
plot(time_opt(1:2:end),host_den(1:2:end),'LineWidth',5,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('\phi 18:3--CBA 18','Interpreter','latex');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')
box on;



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
xline(training_time(5),'k--'); hold on;
plot(time_opt,host_den,'LineWidth',5,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('PSA HS6--PSA H100','Interpreter','latex');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')
box on;


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
xline(training_time(6),'k--'); hold on;
plot(time_opt,host_den,'LineWidth',5,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('PSA HS6--PSA 13-15','Interpreter','latex');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')
box on;



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
xline(training_time(7),'k--'); hold on;
plot(time_opt,host_den,'LineWidth',5,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('PSA HP1--PSA H100','Interpreter','latex');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')
box on;




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
%%
subplot(2,4,8)
xline(training_time(8),'k--'); hold on;
plot(time_opt,host_den,'LineWidth',5,'Color',[0.8,0.8,0.8]);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18);
ylim([1e-3 1e7])
yticks([1e-3 1e-1 1e1 1e3 1e5 1e7]);
xticks(0:2:16);
xlim([0 16]);
title('PSA HP1-PSA 13:15');
ylabel('Host density (genome copies/ml)')
xlabel('Time (hrs)')
box on;