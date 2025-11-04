clc;
clear all;
addpath(genpath(pwd)); % add current directory to path (mcmcstat is included)
addpath('./simulator/');
addpath('./mcmcstat/');


blue1 = [173, 216, 230]/255;
blue2 = [135, 206, 235]/255;
blue3 = [100, 149, 237]/255;
blue4 = [70, 130, 180]/255;
blue5 = [25, 25, 112]/255;
blue6 = [0, 0, 139]/255;


blue = [blue1;blue2;blue3;blue4;blue5;blue6];
Dc_array = [1e4,1e5,1e6,5e6,10e6,20e6];
moi_list = [];
infectious_ratio = [1.72, 2.07, 1.64, 0.01, 3.3, 3.57, 1.87, 1.64]/100;
%% 1

qpcr_host = 1e3*[4.30E+02	2.66E+02	4.71E+02;
2.90E+03	2.27E+03	2.71E+03]*10;


qpcr_virus = 1e3*[1.50E+04	1.50E+04	1.37E+04;
6.01E+08	5.57E+08	5.58E+08]/100;


load('./data_2024/CBA18-2_18_2024.mat');
name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');
format short g


NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(qpcr_host(1,:));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qpcr_virus(1,:));
tvec = 0:0.1:15.8;

moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;

total_virus_seivd_pairwise = cell(1,length(Dc_array));
host_den_seivd_paiwise = cell(1,length(Dc_array));
time_seivd =  cell(1,length(Dc_array));

for i = 1:length(Dc_array)
    [time_seivd_temp, y_series_seivd] = one_step_simulate_seivd(tvec,y0,theta_optimized,NE_optimal,Dc_array(i));
    total_virus_seivd_pairwise_temp = y_series_seivd(end-1,:);
    host_den_seivd_paiwise_temp = sum(y_series_seivd(1:end-2,:));
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end


%subplot(4,4,1)
subplot(2,1,1)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},total_virus_seivd_pairwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('CBA18-2 on CBA 18');
set(gca,'FontSize',20);
xlim([0 16]);
plot([0 945]/60, mean(qpcr_virus'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');


%subplot(4,4,2)
subplot(2,1,2)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},host_den_seivd_paiwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('CBA18-2 on CBA 18');
xlim([0 16]);
plot([0 945]/60, mean(qpcr_host'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');
box on;

%% 2

qpcr_host = 1e3*[1.04E+02	1.14E+02	2.28E+02;
5.37E+02	5.26E+02	5.16E+02];

qpcr_virus = 1e3*[1.11E+04	1.30E+04	1.29E+04;
1.97E+07	1.91E+07	1.90E+07];

moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;

load('./data_2024/CBA18-3_4_2024.mat','theta_optimized');

name = string(labels.phage)+'_'+string(labels.host);
load('parameters.mat','pars');


format short g



data.xdata{1} = time_free_phages;
data.ydata{1} = free_phages;


NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(qpcr_host(1,:));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qpcr_virus(1,:));
tvec = 0:0.1:15.8;


total_virus_seivd_pairwise = cell(1,length(Dc_array));
host_den_seivd_paiwise = cell(1,length(Dc_array));
time_seivd =  cell(1,length(Dc_array));

for i = 1:length(Dc_array)
    [time_seivd_temp, y_series_seivd] = one_step_simulate_seivd(tvec,y0,theta_optimized,NE_optimal,Dc_array(i));
    total_virus_seivd_pairwise_temp = y_series_seivd(end-1,:);
    host_den_seivd_paiwise_temp = sum(y_series_seivd(1:end-2,:));
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end


subplot(4,4,3)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},total_virus_seivd_pairwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('CBA18-3 on CBA 4');
set(gca,'FontSize',20);
xlim([0 16]);
plot([0 945]/60, mean(qpcr_virus'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');


subplot(4,4,4)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},host_den_seivd_paiwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('CBA18-3 on CBA 4');
xlim([0 16]);
plot([0 945]/60, mean(qpcr_host'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');
box on;

%% 3

qpcr_host = 1e3*[3.32E+02	2.26E+02	3.80E+02;
1.45E+02	1.02E+02	7.51E+01];

qpcr_virus = 1e3*[1.08E+04	1.12E+04	1.22E+04;
6.36E+07	5.93E+07	4.98E+07];

moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;

load('./data_2024/CBA18-3_18_2024.mat');


name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');

NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(mean(qpcr_host(1,:)));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qpcr_virus(1,:));

dilution_factor = 100;



tvec = 0:0.01:15.75;

total_virus_seivd_pairwise = cell(1,length(Dc_array));
host_den_seivd_paiwise = cell(1,length(Dc_array));
time_seivd =  cell(1,length(Dc_array));

for i = 1:length(Dc_array)
    [time_seivd_temp, y_series_seivd] = one_step_simulate_seivd(tvec,y0,theta_optimized,NE_optimal,Dc_array(i));
    total_virus_seivd_pairwise_temp = y_series_seivd(end-1,:);
    host_den_seivd_paiwise_temp = sum(y_series_seivd(1:end-2,:));
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end


subplot(4,4,5)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},total_virus_seivd_pairwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('CBA18-3 on CBA 18');
set(gca,'FontSize',20);
xlim([0 16]);
plot([0 945]/60, mean(qpcr_virus'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');


subplot(4,4,6)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},host_den_seivd_paiwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('CBA18-3 on CBA 18');
xlim([0 16]);
plot([0 945]/60, mean(qpcr_host'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');
box on;

%%  4


qpcr_host = 1e3*[1.07E+02	1.15E+02	1.00E+02;
1.28E+02	8.30E+01	2.68E+02];

qpcr_virus = 1e3*[3.86E+05	4.23E+05	4.22E+05;
9.79E+06	9.57E+06	9.60E+06];

moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;

load('./data_2024/CBA38-1_38.mat');

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');

NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(qpcr_host(1,:));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qpcr_virus(1,:));

dilution_factor = 100;


tvec = 0:0.01:15.75;
total_virus_seivd_pairwise = cell(1,length(Dc_array));
host_den_seivd_paiwise = cell(1,length(Dc_array));
time_seivd =  cell(1,length(Dc_array));

for i = 1:length(Dc_array)
    [time_seivd_temp, y_series_seivd] = one_step_simulate_seivd(tvec,y0,theta_optimized,NE_optimal,Dc_array(i));
    total_virus_seivd_pairwise_temp = y_series_seivd(end-1,:);
    host_den_seivd_paiwise_temp = sum(y_series_seivd(1:end-2,:));
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end


subplot(4,4,7)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},total_virus_seivd_pairwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('CBA38-1 on CBA 38');
set(gca,'FontSize',20);
xlim([0 16]);
plot([0 945]/60, mean(qpcr_virus'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');


subplot(4,4,8)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},host_den_seivd_paiwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('CBA38-1 on CBA 38');
xlim([0 16]);
plot([0 945]/60, mean(qpcr_host'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');
box on;

%%  5
qpcr_host = 1e3*[8.79E+01	5.59E+01	2.70E+01;
3.05E+01	3.60E+00	2.01E+01];

qpcr_virus = 1e3*[2.55E+03	3.41E+03	3.51E+03;
8.25E+07	7.98E+07	7.45E+07];

moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;


load('./data_2024/HP1_H100_2024.mat');

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');


format short g


seed = 11;

num_replicates = length(V0_replicates)/3;

NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(mean(qpcr_host(1,:)));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qpcr_virus(1,:));


tvec = 0:0.01:15.75;
total_virus_seivd_pairwise = cell(1,length(Dc_array));
host_den_seivd_paiwise = cell(1,length(Dc_array));
time_seivd =  cell(1,length(Dc_array));

for i = 1:length(Dc_array)
    [time_seivd_temp, y_series_seivd] = one_step_simulate_seivd(tvec,y0,theta_optimized,NE_optimal,Dc_array(i));
    total_virus_seivd_pairwise_temp = y_series_seivd(end-1,:);
    host_den_seivd_paiwise_temp = sum(y_series_seivd(1:end-2,:));
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end


subplot(4,4,9)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},total_virus_seivd_pairwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('HP1 on H100');
set(gca,'FontSize',20);
xlim([0 16]);
plot([0 945]/60, mean(qpcr_virus'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');


subplot(4,4,10)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},host_den_seivd_paiwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('HP1 on H100');
xlim([0 16]);
plot([0 945]/60, mean(qpcr_host'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');
box on;

%% 6


qpcr_host = 1e3*[1.58E+02	9.28E+01	2.59E+02;
8.97E+00	7.95E+00	1.01E+01];

qpcr_virus = 1e3*[4.22E+03	3.60E+03	3.37E+03;
5.54E+07	5.31E+07	5.46E+07];

moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;


load('./data_2024/HP1_13-15_2024.mat');

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');


num_replicates = length(V0_replicates)/3;

NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(mean(qpcr_host(1,:)));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qpcr_virus(1,:));



tvec = 0:0.01:15.75;
total_virus_seivd_pairwise = cell(1,length(Dc_array));
host_den_seivd_paiwise = cell(1,length(Dc_array));
time_seivd =  cell(1,length(Dc_array));

for i = 1:length(Dc_array)
    [time_seivd_temp, y_series_seivd] = one_step_simulate_seivd(tvec,y0,theta_optimized,NE_optimal,Dc_array(i));
    total_virus_seivd_pairwise_temp = y_series_seivd(end-1,:);
    host_den_seivd_paiwise_temp = sum(y_series_seivd(1:end-2,:));
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end


subplot(4,4,11)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},total_virus_seivd_pairwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('HP1 on PSA 13-15');
set(gca,'FontSize',20);
xlim([0 16]);
plot([0 945]/60, mean(qpcr_virus'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');


subplot(4,4,12)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},host_den_seivd_paiwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('HP1 on PSA 13-15');
xlim([0 16]);
plot([0 945]/60, mean(qpcr_host'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');
box on;

%%  7

qpcr_host = 1e3*[1.16E+01	6.26E+01	6.66E+01;
9.89E+02	2.36E+02	9.40E+01];

qpcr_virus = 1e3*[1.05E+04	6.57E+06	1.09E+04;
2.08E+08	1.68E+08	1.85E+08];


moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;


load('./data_2024/HS6_H100_2024.mat');

name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');



num_replicates = length(V0_replicates)/3;


NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(mean(qpcr_host(1,:)));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qpcr_virus(1,:));


tvec = 0:0.01:15.75;
total_virus_seivd_pairwise = cell(1,length(Dc_array));
host_den_seivd_paiwise = cell(1,length(Dc_array));
time_seivd =  cell(1,length(Dc_array));

for i = 1:length(Dc_array)
    [time_seivd_temp, y_series_seivd] = one_step_simulate_seivd(tvec,y0,theta_optimized,NE_optimal,Dc_array(i));
    total_virus_seivd_pairwise_temp = y_series_seivd(end-1,:);
    host_den_seivd_paiwise_temp = sum(y_series_seivd(1:end-2,:));
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end


subplot(4,4,13)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},total_virus_seivd_pairwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('HS6 on H100');
set(gca,'FontSize',20);
xlim([0 16]);
plot([0 945]/60, mean(qpcr_virus'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');


subplot(4,4,14)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},host_den_seivd_paiwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('HS6 on H100');
xlim([0 16]);
plot([0 945]/60, mean(qpcr_host'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');
box on;

%%  8


qpcr_host = 1e3*[9.02E+01	4.22E+02	1.03E+02;
6.25E+02	7.06E+02	1.30E+03];

qpcr_virus = 1e3*[1.02E+04	1.07E+04	1.21E+04;
1.03E+08	9.02E+07	9.36E+07];


moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;


load('./data_2024/HS6_13-15_2024.mat');


name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');


format short g



num_replicates = length(V0_replicates)/3;




NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(mean(qpcr_host(1,:)));
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(qpcr_virus(1,:));


tvec = 0:0.01:15.75;
total_virus_seivd_pairwise = cell(1,length(Dc_array));
host_den_seivd_paiwise = cell(1,length(Dc_array));
time_seivd =  cell(1,length(Dc_array));

for i = 1:length(Dc_array)
    [time_seivd_temp, y_series_seivd] = one_step_simulate_seivd(tvec,y0,theta_optimized,NE_optimal,Dc_array(i));
    total_virus_seivd_pairwise_temp = y_series_seivd(end-1,:);
    host_den_seivd_paiwise_temp = sum(y_series_seivd(1:end-2,:));
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end


subplot(4,4,15)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},total_virus_seivd_pairwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)')
ylabel('Virus density (/ml)')
title('HS6 on 13-15');
set(gca,'FontSize',20);
xlim([0 16]);
plot([0 945]/60, mean(qpcr_virus'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');


subplot(4,4,16)
for i = 1:length(Dc_array)
hold on;
plot(time_seivd{i},host_den_seivd_paiwise{i},'LineWidth',2,'Color',blue(i,:));
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('HS6 on 13-15');
xlim([0 16]);
plot([0 945]/60, mean(qpcr_host'), 'Marker','square','MarkerEdgeColor','k','MarkerFaceColor','k','LineStyle','none');
box on;


% Display the legend for this subplot
legend('D_c = 1e4 cells/ml', ...
    'D_c = 1e5 cells/ml', ...
     'D_c = 1e6 cells/ml', ...
      'D_c = 5e6 cells/ml', ...
       'D_c = 10e6 cells/ml', ...
        'D_c = 20e6 cells/ml')


legend('show');

%%

moi_list

moi_list.*infectious_ratio