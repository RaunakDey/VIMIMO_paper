clc;
clear all;
addpath(genpath(pwd)); % add current directory to path (mcmcstat is included)
addpath('./simulator/');
addpath('./mcmcstat/');

%%




%% 1

qprc_host = 1e3*[4.30E+02	2.66E+02	4.71E+02;
2.90E+03	2.27E+03	2.71E+03];


qprc_virus = 1e3*[1.50E+04	1.50E+04	1.37E+04;
6.01E+08	5.57E+08	5.58E+08];




load('./data_2024/CBA18-2_18_2024.mat');

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');


format short g

num_replicates = 1;


NE_optimal = round(theta_optimized(5));
clear y0

%S0_replicates = qprc_host(1,:);
%V0_replicates = qprc_virus(1,:);

%S0_replicates = 2.4e6;
%V0_replicates = 2.4e5;

y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);

dilution_factor = 100;

tvec = 0:0.1:15;
[time_opt,y_series_opt] = one_step_simulate(tvec,y0,theta_optimized,NE_optimal,dilution_factor);

%%%%% new models
prob_effective_infection = 0.5;
prob_lysis_reset = 0.1;

theta_optimized(6) = prob_effective_infection; %prob effective inhibition
theta_optimized(7) = prob_lysis_reset; %prob lysis reset
NE_optimal = round(theta_optimized(5));
clear y0

%S0_replicates = qprc_host(1,:);
%V0_replicates = qprc_virus(1,:);

%S0_replicates = 2.4e6;
%V0_replicates = 2.4e5;

y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);
dilution_factor = 100;

[time_ineffi_adsorb,y_series_ineffi_adsorb] = one_step_simulate_infec_vary_adsorbs_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_adsorption_model = y_series_ineffi_adsorb(end,:)+y_series_ineffi_adsorb(end-1,:);

[time_ineffi_infection, y_series_ineffi_infection] = one_step_simulate_infec_vary_infec_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_infection_model = y_series_ineffi_infection(end,:)+y_series_ineffi_infection(end-1,:);

[time_lysis_res, y_series_lysis_res] = one_step_simulate_lysis_reset(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_lysis_res = y_series_lysis_res(end,:);



subplot(4,4,1)

plot(time_opt,y_series_opt(end,:),'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,total_virus_ineffi_adsorption_model,'-r','LineWidth',2);
plot(time_ineffi_infection,total_virus_ineffi_infection_model,'-g','LineWidth',2);
plot(time_lysis_res,total_virus_lysis_res,'-b','LineWidth',2);

set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('CBA18-2 on CBA 18');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'mo'); hold on;
    end
end
xlim([0 15.75]);
%ylim([1e3 1e9]);
yticks([1e3 1e5 1e7 1e9]);
plot([0 945]/60, mean(qprc_virus'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');



host_den = sum(y_series_opt(1:end-1,:));
host_den_ineff_adsorb = sum(y_series_ineffi_adsorb(1:end-2,:));
host_den_ineff_infection =  sum(y_series_ineffi_infection(1:end-2,:));
host_den_lysis_res = sum(y_series_lysis_res(1:end-1,:));

subplot(4,4,2)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,host_den_ineff_adsorb,'-r','LineWidth',2);hold on;
plot(time_ineffi_infection,host_den_ineff_infection,'-g','LineWidth',2);hold on;
plot(time_lysis_res,host_den_lysis_res,'-b','LineWidth',2);hold on;

set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('CBA18-2 on CBA 18');
xlim([0 15.75]);
%ylim([1e-1 1e7]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);



plot([0 945]/60, mean(qprc_host'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');


%% 2


load('./data_2024/CBA18-3_4_2024.mat');

name = string(labels.phage)+'_'+string(labels.host);
load('parameters.mat','pars');


format short g



data.xdata{1} = time_free_phages;
data.ydata{1} = free_phages;


NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);

dilution_factor = 100;

tvec = 0:0.01:15.75;
[time_opt,y_series_opt] = one_step_simulate(tvec,y0,theta_optimized,NE_optimal,dilution_factor);

%%%%% new models
prob_effective_infection = 0.5;
prob_lysis_reset = 0.1;

theta_optimized(6) = prob_effective_infection; %prob effective inhibition
theta_optimized(7) = prob_lysis_reset; %prob lysis reset
NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);
dilution_factor = 100;

[time_ineffi_adsorb,y_series_ineffi_adsorb] = one_step_simulate_infec_vary_adsorbs_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_adsorption_model = y_series_ineffi_adsorb(end,:)+y_series_ineffi_adsorb(end-1,:);

[time_ineffi_infection, y_series_ineffi_infection] = one_step_simulate_infec_vary_infec_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_infection_model = y_series_ineffi_infection(end,:)+y_series_ineffi_infection(end-1,:);


[time_lysis_res, y_series_lysis_res] = one_step_simulate_lysis_reset(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_lysis_res = y_series_lysis_res(end,:);


host_den = sum(y_series_opt(1:end-1,:));
host_den_ineff_adsorb = sum(y_series_ineffi_adsorb(1:end-2,:));
host_den_ineff_infection =  sum(y_series_ineffi_infection(1:end-2,:));
host_den_lysis_res = sum(y_series_lysis_res(1:end-1,:));

subplot(4,4,3)

plot(time_opt,y_series_opt(end,:),'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,total_virus_ineffi_adsorption_model,'-r','LineWidth',2);
plot(time_ineffi_infection,total_virus_ineffi_infection_model,'-g','LineWidth',2);
plot(time_lysis_res,total_virus_lysis_res,'-b','LineWidth',2);

set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('CBA 18:3 on CBA 4');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'mo'); hold on;
    end
end
xlim([0 15.75]);
ylim([1e3 1e9]);
yticks([1e3 1e5 1e7 1e9])




subplot(4,4,4)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,host_den_ineff_adsorb,'-r','LineWidth',2);hold on;
plot(time_ineffi_infection,host_den_ineff_infection,'-g','LineWidth',2);hold on;
plot(time_lysis_res,host_den_lysis_res,'-b','LineWidth',2);hold on;

set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('CBA 18:3 on CBA 4');
xlim([0 15.75]);
ylim([1e-1 1e7]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);


%% 3


load('./data_2024/CBA18-3_18_2024.mat');


name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');

NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);

dilution_factor = 100;



tvec = 0:0.01:15.75;
[time_opt,y_series_opt] = one_step_simulate(tvec,y0,theta_optimized,NE_optimal,dilution_factor);

%%%%% new models
prob_effective_infection = 0.5;
prob_lysis_reset = 0.1;

theta_optimized(6) = prob_effective_infection; %prob effective inhibition
theta_optimized(7) = prob_lysis_reset; %prob lysis reset
NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);
dilution_factor = 100;

[time_ineffi_adsorb,y_series_ineffi_adsorb] = one_step_simulate_infec_vary_adsorbs_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_adsorption_model = y_series_ineffi_adsorb(end,:)+y_series_ineffi_adsorb(end-1,:);

[time_ineffi_infection, y_series_ineffi_infection] = one_step_simulate_infec_vary_infec_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_infection_model = y_series_ineffi_infection(end,:)+y_series_ineffi_infection(end-1,:);


[time_lysis_res, y_series_lysis_res] = one_step_simulate_lysis_reset(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_lysis_res = y_series_lysis_res(end,:);


host_den = sum(y_series_opt(1:end-1,:));
host_den_ineff_adsorb = sum(y_series_ineffi_adsorb(1:end-2,:));
host_den_ineff_infection =  sum(y_series_ineffi_infection(1:end-2,:));
host_den_lysis_res = sum(y_series_lysis_res(1:end-1,:));

subplot(4,4,5)


plot(time_opt,y_series_opt(end,:),'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,total_virus_ineffi_adsorption_model,'-r','LineWidth',2);
plot(time_ineffi_infection,total_virus_ineffi_infection_model,'-g','LineWidth',2);
plot(time_lysis_res,total_virus_lysis_res,'-b','LineWidth',2);

set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('CBA18:3 on CBA18');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'mo'); hold on;
    end
end
xlim([0 15.75]);
ylim([1e3 1e9]);
yticks([1e3 1e5 1e7 1e9])


subplot(4,4,6)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,host_den_ineff_adsorb,'-r','LineWidth',2);hold on;
plot(time_ineffi_infection,host_den_ineff_infection,'-g','LineWidth',2);hold on;
plot(time_lysis_res,host_den_lysis_res,'-b','LineWidth',2);hold on;

set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('CBA18:3 on CBA18');
xlim([0 15.75]);
ylim([1e-1 1e7]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);


%%  4


load('./data_2024/CBA38-1_38.mat');

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');

NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);

dilution_factor = 100;


tvec = 0:0.01:15.75;
[time_opt,y_series_opt] = one_step_simulate(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
%%%%% new models
prob_effective_infection = 0.5;
prob_lysis_reset = 0.1;

theta_optimized(6) = prob_effective_infection; %prob effective inhibition
theta_optimized(7) = prob_lysis_reset; %prob lysis reset
NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);
dilution_factor = 100;

[time_ineffi_adsorb,y_series_ineffi_adsorb] = one_step_simulate_infec_vary_adsorbs_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_adsorption_model = y_series_ineffi_adsorb(end,:)+y_series_ineffi_adsorb(end-1,:);

[time_ineffi_infection, y_series_ineffi_infection] = one_step_simulate_infec_vary_infec_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_infection_model = y_series_ineffi_infection(end,:)+y_series_ineffi_infection(end-1,:);


[time_lysis_res, y_series_lysis_res] = one_step_simulate_lysis_reset(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_lysis_res = y_series_lysis_res(end,:);


host_den = sum(y_series_opt(1:end-1,:));
host_den_ineff_adsorb = sum(y_series_ineffi_adsorb(1:end-2,:));
host_den_ineff_infection =  sum(y_series_ineffi_infection(1:end-2,:));
host_den_lysis_res = sum(y_series_lysis_res(1:end-1,:));


subplot(4,4,7)
plot(time_opt,y_series_opt(end,:),'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,total_virus_ineffi_adsorption_model,'-r','LineWidth',2);
plot(time_ineffi_infection,total_virus_ineffi_infection_model,'-g','LineWidth',2);
plot(time_lysis_res,total_virus_lysis_res,'-b','LineWidth',2);

set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('CBA38:1 on CBA38');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'mo'); hold on;
    end
end
xlim([0 15.75]);
ylim([1e3 1e9]);
yticks([1e3 1e5 1e7 1e9]);


subplot(4,4,8)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,host_den_ineff_adsorb,'-r','LineWidth',2);hold on;
plot(time_ineffi_infection,host_den_ineff_infection,'-g','LineWidth',2);hold on;
plot(time_lysis_res,host_den_lysis_res,'-b','LineWidth',2);hold on;

set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('CBA38:1 on CBA38');
xlim([0 15.75]);
ylim([1e-1 1e7]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);


%%  5


load('./data_2024/HP1_H100_2024.mat');

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');


format short g


seed = 11;

num_replicates = length(V0_replicates)/3;

NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);

dilution_factor = 100;
tvec = 0:0.01:15.75;
[time_opt,y_series_opt] = one_step_simulate(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
%%%%% new models
prob_effective_infection = 0.25;
prob_lysis_reset = 0.1;

theta_optimized(6) = prob_effective_infection; %prob effective inhibition
theta_optimized(7) = prob_lysis_reset; %prob lysis reset
NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);
dilution_factor = 100;

[time_ineffi_adsorb,y_series_ineffi_adsorb] = one_step_simulate_infec_vary_adsorbs_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_adsorption_model = y_series_ineffi_adsorb(end,:)+y_series_ineffi_adsorb(end-1,:);

[time_ineffi_infection, y_series_ineffi_infection] = one_step_simulate_infec_vary_infec_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_infection_model = y_series_ineffi_infection(end,:)+y_series_ineffi_infection(end-1,:);


[time_lysis_res, y_series_lysis_res] = one_step_simulate_lysis_reset(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_lysis_res = y_series_lysis_res(end,:);


host_den = sum(y_series_opt(1:end-1,:));
host_den_ineff_adsorb = sum(y_series_ineffi_adsorb(1:end-2,:));
host_den_ineff_infection =  sum(y_series_ineffi_infection(1:end-2,:));
host_den_lysis_res = sum(y_series_lysis_res(1:end-1,:));



subplot(4,4,9)
plot(time_opt,y_series_opt(end,:),'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,total_virus_ineffi_adsorption_model,'-r','LineWidth',2);
plot(time_ineffi_infection,total_virus_ineffi_infection_model,'-g','LineWidth',2);
plot(time_lysis_res,total_virus_lysis_res,'-b','LineWidth',2);


set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('HP1 on H100');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'mo'); hold on;
    end
end
xlim([0 15.75]);
ylim([1e3 1e9]);
yticks([1e3 1e5 1e7 1e9]);




subplot(4,4,10)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,host_den_ineff_adsorb,'-r','LineWidth',2);hold on;
plot(time_ineffi_infection,host_den_ineff_infection,'-g','LineWidth',2);hold on;
plot(time_lysis_res,host_den_lysis_res,'-b','LineWidth',2);hold on;

set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('HP1 on H100');
xlim([0 15.75]);
ylim([1e-1 1e7]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);


%% 6

load('./data_2024/HP1_13-15_2024.mat');

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');


num_replicates = length(V0_replicates)/3;

NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);

dilution_factor = 100;

tvec = 0:0.01:15.75;
[time_opt,y_series_opt] = one_step_simulate(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
%%%%% new models
prob_effective_infection = 0.25;
prob_lysis_reset = 0.1;

theta_optimized(6) = prob_effective_infection; %prob effective inhibition
theta_optimized(7) = prob_lysis_reset; %prob lysis reset
NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);
dilution_factor = 100;

[time_ineffi_adsorb,y_series_ineffi_adsorb] = one_step_simulate_infec_vary_adsorbs_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_adsorption_model = y_series_ineffi_adsorb(end,:)+y_series_ineffi_adsorb(end-1,:);

[time_ineffi_infection, y_series_ineffi_infection] = one_step_simulate_infec_vary_infec_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_infection_model = y_series_ineffi_infection(end,:)+y_series_ineffi_infection(end-1,:);


[time_lysis_res, y_series_lysis_res] = one_step_simulate_lysis_reset(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_lysis_res = y_series_lysis_res(end,:);


host_den = sum(y_series_opt(1:end-1,:));
host_den_ineff_adsorb = sum(y_series_ineffi_adsorb(1:end-2,:));
host_den_ineff_infection =  sum(y_series_ineffi_infection(1:end-2,:));
host_den_lysis_res = sum(y_series_lysis_res(1:end-1,:));


subplot(4,4,11)
plot(time_opt,y_series_opt(end,:),'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,total_virus_ineffi_adsorption_model,'-r','LineWidth',2);
plot(time_ineffi_infection,total_virus_ineffi_infection_model,'-g','LineWidth',2);
plot(time_lysis_res,total_virus_lysis_res,'-b','LineWidth',2);

set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('HP1 on 13-15');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'mo'); hold on;
    end
end
xlim([0 15.75]);
ylim([1e3 1e9]);
yticks([1e3 1e5 1e7 1e9]);



subplot(4,4,12)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,host_den_ineff_adsorb,'-r','LineWidth',2);hold on;
plot(time_ineffi_infection,host_den_ineff_infection,'-g','LineWidth',2);hold on;
plot(time_lysis_res,host_den_lysis_res,'-b','LineWidth',2);hold on;

set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('HP1 on 13-15');
xlim([0 15.75]);
ylim([1e-1 1e7]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);


%%  7


load('./data_2024/HS6_H100_2024.mat');

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');



num_replicates = length(V0_replicates)/3;


NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);

dilution_factor = 100;

tvec = 0:0.01:15.75;
[time_opt,y_series_opt] = one_step_simulate(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
%%%%% new models
prob_effective_infection = 0.15;
prob_lysis_reset = 0.1;

theta_optimized(6) = prob_effective_infection; %prob effective inhibition
theta_optimized(7) = prob_lysis_reset; %prob lysis reset
NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);
dilution_factor = 100;

[time_ineffi_adsorb,y_series_ineffi_adsorb] = one_step_simulate_infec_vary_adsorbs_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_adsorption_model = y_series_ineffi_adsorb(end,:)+y_series_ineffi_adsorb(end-1,:);

[time_ineffi_infection, y_series_ineffi_infection] = one_step_simulate_infec_vary_infec_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_infection_model = y_series_ineffi_infection(end,:)+y_series_ineffi_infection(end-1,:);


[time_lysis_res, y_series_lysis_res] = one_step_simulate_lysis_reset(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_lysis_res = y_series_lysis_res(end,:);


host_den = sum(y_series_opt(1:end-1,:));
host_den_ineff_adsorb = sum(y_series_ineffi_adsorb(1:end-2,:));
host_den_ineff_infection =  sum(y_series_ineffi_infection(1:end-2,:));
host_den_lysis_res = sum(y_series_lysis_res(1:end-1,:));



subplot(4,4,13)
plot(time_opt,y_series_opt(end,:),'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,total_virus_ineffi_adsorption_model,'-r','LineWidth',2);
plot(time_ineffi_infection,total_virus_ineffi_infection_model,'-g','LineWidth',2);
plot(time_lysis_res,total_virus_lysis_res,'-b','LineWidth',2);

set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('HS6 on H100');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'mo'); hold on;
    end
end
xlim([0 15.75]);
ylim([1e3 1e9]);
yticks([1e3 1e5 1e7 1e9]);



subplot(4,4,14)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,host_den_ineff_adsorb,'-r','LineWidth',2);hold on;
plot(time_ineffi_infection,host_den_ineff_infection,'-g','LineWidth',2);hold on;
plot(time_lysis_res,host_den_lysis_res,'-b','LineWidth',2);hold on;

set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('HS6 on H100');
xlim([0 15.75]);
ylim([1e-1 1e7]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);



%%  8


load('./data_2024/HS6_13-15_2024.mat');


name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');


format short g



num_replicates = length(V0_replicates)/3;




NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);

dilution_factor = 100;

tvec = 0:0.01:15.75;
[time_opt,y_series_opt] = one_step_simulate(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
%%%%% new models
prob_effective_infection = 0.15;
prob_lysis_reset = 0.1;

theta_optimized(6) = prob_effective_infection; %prob effective inhibition
theta_optimized(7) = prob_lysis_reset; %prob lysis reset
NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);
dilution_factor = 100;

[time_ineffi_adsorb,y_series_ineffi_adsorb] = one_step_simulate_infec_vary_adsorbs_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_adsorption_model = y_series_ineffi_adsorb(end,:)+y_series_ineffi_adsorb(end-1,:);

[time_ineffi_infection, y_series_ineffi_infection] = one_step_simulate_infec_vary_infec_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_ineffi_infection_model = y_series_ineffi_infection(end,:)+y_series_ineffi_infection(end-1,:);


[time_lysis_res, y_series_lysis_res] = one_step_simulate_lysis_reset(tvec,y0,theta_optimized,NE_optimal,dilution_factor);
total_virus_lysis_res = y_series_lysis_res(end,:);


host_den = sum(y_series_opt(1:end-1,:));
host_den_ineff_adsorb = sum(y_series_ineffi_adsorb(1:end-2,:));
host_den_ineff_infection =  sum(y_series_ineffi_infection(1:end-2,:));
host_den_lysis_res = sum(y_series_lysis_res(1:end-1,:));


subplot(4,4,15)

% Plot the virus density
plot(time_opt,y_series_opt(end,:),'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,total_virus_ineffi_adsorption_model,'-r','LineWidth',2);
plot(time_ineffi_infection,total_virus_ineffi_infection_model,'-g','LineWidth',2);
plot(time_lysis_res,total_virus_lysis_res,'-b','LineWidth',2);

hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Virus density (/ml)');
title('HS6 on 13-15');

set(gca,'FontSize',20);

% Plot the virus data points, but without adding them to the legend
for i = 1:num_replicates
    for j = 1:3
        plot(data.xdata{i}./60, data.ydata{i}(:,j), 'mo');
    end
end

% Add one of the data series to the legend
plot(data.xdata{i}./60, data.ydata{i}(:,j), 'mo', 'DisplayName', 'Virus Data');

xlim([0 15.75]);
ylim([1e3 1e9]);
yticks([1e3 1e5 1e7 1e9]);

% Display the legend for this subplot
legend('Baseline model','Ineffective viral adsorption','Ineffective viral infection','Lysis reset model','Virus Data','AutoUpdate','off');
legend('show');

% Plot the host density in the second subplot
host_den = sum(y_series_opt(1:end-1,:));
subplot(4,4,16);
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
plot(time_ineffi_adsorb,host_den_ineff_adsorb,'-r','LineWidth',2);hold on;
plot(time_ineffi_infection,host_den_ineff_infection,'-g','LineWidth',2);hold on;
plot(time_lysis_res,host_den_lysis_res,'-b','LineWidth',2);hold on;

hold on;
set(gca, 'YScale', 'log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca, 'FontSize', 20);
title('HS6 on 13-15');
xlim([0 15.75]);
ylim([1e-1 1e7]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);

% Display the legend for this subplot
%legend('Virus Density Baseline','Ineffective viral adsorption','Ineffective viral infection','Virus Data','AutoUpdate','off');

%legend('show');
