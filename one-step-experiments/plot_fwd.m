clc;
clear all;
addpath(genpath(pwd)); % add current directory to path (mcmcstat is included)
addpath('./simulator/');
addpath('./mcmcstat/');




%% 

load('./data_2024/CBA18-2_18_2024.mat');

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');


format short g

num_replicates = 1;


NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);

dilution_factor = 100;

tvec = 0:0.01:15.75;
[time_opt,y_series_opt] = one_step_simulate(tvec,y0,theta_optimized,NE_optimal,dilution_factor);

figure
subplot(2,1,1)

plot(time_opt,y_series_opt(end,:),'-g','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('CBA18-2 on CBA 18');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'bo'); hold on;
    end
end




host_den = sum(y_series_opt(1:end-1,:));
subplot(2,1,2)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);



%% 


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

figure
subplot(2,1,1)

plot(time_opt,y_series_opt(end,:),'-g','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('CBA 18:3 on CBA 4');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'bo'); hold on;
    end
end



host_den = sum(y_series_opt(1:end-1,:));
subplot(2,1,2)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);

%%


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

figure
subplot(2,1,1)

plot(time_opt,y_series_opt(end,:),'-g','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('CBA18:3 on CBA18');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'bo'); hold on;
    end
end



host_den = sum(y_series_opt(1:end-1,:));
subplot(2,1,2)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);


%%  


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


figure
subplot(2,1,1)

plot(time_opt,y_series_opt(end,:),'-g','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('CBA38:1 on CBA38');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'bo'); hold on;
    end
end



host_den = sum(y_series_opt(1:end-1,:));
subplot(2,1,2)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);



%%  


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


figure
subplot(2,1,1)

plot(time_opt,y_series_opt(end,:),'-g','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('HP1 on H100');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'bo'); hold on;
    end
end



host_den = sum(y_series_opt(1:end-1,:));
subplot(2,1,2)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);


%% 

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


figure
subplot(2,1,1)

plot(time_opt,y_series_opt(end,:),'-g','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('HP1 on 13-15');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'bo'); hold on;
    end
end



host_den = sum(y_series_opt(1:end-1,:));
subplot(2,1,2)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);


%% 


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


figure
subplot(2,1,1)

plot(time_opt,y_series_opt(end,:),'-g','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('HS6 on H100');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'bo'); hold on;
    end
end



host_den = sum(y_series_opt(1:end-1,:));
subplot(2,1,2)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);



%% 


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


figure
subplot(2,1,1)

plot(time_opt,y_series_opt(end,:),'-g','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('HS6 on H100');
set(gca,'FontSize',20);
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'bo'); hold on;
    end
end



host_den = sum(y_series_opt(1:end-1,:));
subplot(2,1,2)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);



