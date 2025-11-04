clc;
clear all;
addpath(genpath(pwd)); % add current directory to path (mcmcstat is included)
addpath('./../simulator/');
addpath('./../mcmcstat/');

gray1 = [220 220 220]/255;   % light gray
gray2 = [169 169 169]/255;   % dark gray
gray3 = [128 128 128]/255;   % medium gray
%gray4 = [70 70 70]/255;      % very dark gray
gray5 = [0 0 0]/255;         % black

gray = [gray1; gray2; gray3; gray5];
Dc_array = [1e5,1e6,5e6,1e15];
%Dc_array = [1e5];
moi_list = [];
infectious_ratio = [1.72, 2.07, 1.64, 0.01, 3.3, 3.57, 1.87, 1.64]/100;

cutoff_den = 0.99;

%%

set(0,'DefaultFigureRenderer','painters');
fig = figure('Renderer','painters','Color','w');
set(fig,'GraphicsSmoothing','off');


%% no phage control


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


%% 1) 18:2 -- 18


qpcr_host = 1e3*[5.46E+03	6.68E+03	3.78E+03;
9.15E+03	9.78E+03	7.62E+03;
4.01E+04	3.70E+04	4.32E+04];

qpcr_virus = 1e3*[1.73E+03	1.89E+03	2.53E+03;
1.56E+04	1.10E+04	3.18E+04;
1.16E+08	1.26E+08	1.14E+08];


load('./../data_2024/CBA18-2_18_2024.mat');
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
    host_den_seivd_paiwise_temp(host_den_seivd_paiwise_temp < cutoff_den) = 0;
    
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end



subplot(4,4,2)

for i = 1:length(Dc_array)
    hold on;
    if i == length(Dc_array)
        plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
    else
        plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:));
    end
end

hold on;
set(gca,'YScale','log');
xlabel('Time (hrs)')
ylabel('Phage/ml')
%title('\phi18:2 – CBA 18');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
xlim([0 16]);

mean_virus = mean(qpcr_virus,2);             % mean across replicates (rows are time points)
err = std(qpcr_virus,0,2);          % std across replicates
errorbar([0 175 945]/60, mean_virus, err, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'LineStyle','none', ...
    'Color','r','LineWidth',1); 


ylim([1e4 5e11]);
yticks([1e4, 1e6, 1e8, 1e10]);
xticks([0,2,4,6,8,10,12,14,16]);
box on; xticks([0 2 4 6 8 10 12 14 16])


subplot(4,4,1)


for i = 1:length(Dc_array)
    hold on;
    if i == length(Dc_array)
     

        plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
    else
        plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:));
    end
end



hold on;
set(gca,'YScale','log');
xlabel('Time (hrs)');
ylabel('Host/ml');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
%title('\phi18:2 – CBA 18');
xlim([0 16]);
mean_host = mean(qpcr_host,2);             % mean across replicates (rows are time points)
err = std(qpcr_host,0,2);          % std across replicates
errorbar([0 175 945]/60, mean_host, err, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'LineStyle','none', ...
    'Color','b','LineWidth',1); 

mean_host_nophage = mean(cba18_nophage,2); 
std_host_nophage = std(cba18_nophage,0,2);

errorbar([0 175 945]/60, mean_host_nophage, std_host_nophage, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'b', ...
    'LineStyle','none', ...
    'Color','b','LineWidth',1); 


xticks([0,2,4,6,8,10,12,14,16]);
yticks([1e-6 1e-3 1e0 1e3 1e6 1e9]);
box on; xticks([0 2 4 6 8 10 12 14 16])

%% 2) 18:3 – 4

qpcr_host = 1e3*[3.23E+03	4.18E+03	2.56E+03;
2.20E+03	2.83E+03	1.72E+03;
9.39E+03	3.41E+03	5.98E+03];

qpcr_virus = 1e3*[4.46E+02	4.37E+02	4.40E+02;
5.33E+02	5.30E+02	5.15E+02;
3.69E+05	3.89E+05	4.54E+05];

moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;

load('./../data_2024/CBA18-3_4_2024.mat');

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
    host_den_seivd_paiwise_temp(host_den_seivd_paiwise_temp < cutoff_den) = 0;
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end

subplot(4,4,4)
for i = 1:length(Dc_array)
    hold on;
    if i == length(Dc_array)
        plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
    else
        plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:));
    end
end
hold on;
set(gca,'YScale','log');
xlabel('Time (hrs)')
ylabel('Phage/ml')
%title('\phi18:3 – CBA 4');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
xlim([0 16]);

mean_virus = mean(qpcr_virus,2);             % mean across replicates (rows are time points)
err = std(qpcr_virus,0,2);          % std across replicates
errorbar([0 175 945]/60, mean_virus, err, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'LineStyle','none', ...
    'Color','r','LineWidth',1); 

ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]); xticks([0,2,4,6,8,10,12,14,16]);
box on; xticks([0 2 4 6 8 10 12 14 16])


subplot(4,4,3)

for i = 1:length(Dc_array)
    hold on;
    if i == length(Dc_array)
     

        plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
    else
        plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:));
    end
end



hold on;
set(gca,'YScale','log');
xlabel('Time (hrs)');
ylabel('Host/ml');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
%title('\phi18:3 – CBA 4');
xlim([0 16]);
mean_host = mean(qpcr_host,2);             % mean across replicates (rows are time points)
err = std(qpcr_host,0,2);          % std across replicates
errorbar([0 175 945]/60, mean_host, err, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'LineStyle','None',...
    'Color','b','LineWidth',1); 

mean_host_nophage = mean(cba4_nophage,2); 
std_host_nophage = std(cba4_nophage,0,2);

errorbar([0 175 945]/60, mean_host_nophage, std_host_nophage, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'b', ...
    'LineStyle','none', ...
    'Color','b','LineWidth',1); 



box on; xticks([0 2 4 6 8 10 12 14 16])

%% 3) 18:3 – 18

qpcr_host = 1e3*[7.87E+03	8.00E+03	6.09E+03;
1.04E+04	1.08E+04	6.60E+03;
1.16E+03	1.36E+03	1.58E+03];

qpcr_virus = 1e3*[4.54E+02	4.71E+02	5.45E+02;
3.78E+03	2.70E+03	2.99E+03;
1.57E+06	1.66E+06	1.37E+06];


moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;

load('./../data_2024/CBA18-3_18_2024.mat');


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
        host_den_seivd_paiwise_temp(host_den_seivd_paiwise_temp < cutoff_den) = 0;
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end


subplot(4,4,6)
for i = 1:length(Dc_array)
    hold on;
    if i == length(Dc_array)
        plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
    else
        plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:));
    end
end


set(gca,'YScale','log');
xlabel('Time (hrs)')
ylabel('Phage/ml')
%title('\phi18:3 – CBA 18');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
xlim([0 16]);

mean_virus = mean(qpcr_virus,2);             % mean across replicates (rows are time points)
err = std(qpcr_virus,0,2);          % std across replicates
errorbar([0 175 945]/60, mean_virus, err, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'LineStyle','none', ...
    'Color','r','LineWidth',1); 


ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]); xticks([0,2,4,6,8,10,12,14,16]);
box on; xticks([0 2 4 6 8 10 12 14 16])

subplot(4,4,5)
for i = 1:length(Dc_array)
    hold on;
    if i == length(Dc_array)
     

        plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
    else
        plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:));
    end
end



hold on;
set(gca,'YScale','log');
xlabel('Time (hrs)');
ylabel('Host/ml');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
%title('\phi18:3 – CBA 18');
xlim([0 16]);
mean_host = mean(qpcr_host,2);             % mean across replicates (rows are time points)
err = std(qpcr_host,0,2);          % std across replicates
errorbar([0 175 945]/60, mean_host, err, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'LineStyle','None',...
    'Color','b','LineWidth',1); 

mean_host_nophage = mean(cba18_nophage,2); 
std_host_nophage = std(cba18_nophage,0,2);

errorbar([0 175 945]/60, mean_host_nophage, std_host_nophage, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'b', ...
    'LineStyle','none', ...
    'Color','b','LineWidth',1); 



box on; xticks([0 2 4 6 8 10 12 14 16])
yticks([1e-6 1e-3 1e0 1e3 1e6 1e9]);

%%  4) 38:1 – 38


qpcr_host = 1e3*[9.01E+03	8.17E+03	8.13E+03;
1.09E+04	1.35E+04	9.19E+03;
1.18E+03	1.42E+03	1.86E+03];

qpcr_virus = 1e3*[2.85E+02	2.73E+02	2.90E+02;
2.65E+03	2.48E+03	2.23E+03;
2.81E+05	3.08E+05	3.20E+05];

moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;

load('./../data_2024/CBA38-1_38.mat');

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
        host_den_seivd_paiwise_temp(host_den_seivd_paiwise_temp < cutoff_den) = 0;
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end


subplot(4,4,8)
for i = 1:length(Dc_array)
    hold on;
    if i == length(Dc_array)
        plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
    else
        plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:));
    end
end

hold on;
set(gca,'YScale','log');
xlabel('Time (hrs)')
ylabel('Phage/ml')
%title('\phi38:1 – CBA 38');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
xlim([0 16]);

mean_virus = mean(qpcr_virus,2);             % mean across replicates (rows are time points)
err = std(qpcr_virus,0,2);          % std across replicates
errorbar([0 175 945]/60, mean_virus, err, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'LineStyle','none', ...
    'Color','r','LineWidth',1); 
ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]); xticks([0,2,4,6,8,10,12,14,16]);
box on; xticks([0 2 4 6 8 10 12 14 16])


subplot(4,4,7)
for i = 1:length(Dc_array)
    hold on;
    if i == length(Dc_array)
     

        plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
    else
        plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:));
    end
end



hold on;
set(gca,'YScale','log');
xlabel('Time (hrs)');
ylabel('Host/ml');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
%title('\phi38:1 – CBA 38');
xlim([0 16]);

mean_host = mean(qpcr_host,2);             % mean across replicates (rows are time points)
err = std(qpcr_host,0,2);          % std across replicates
errorbar([0 175 945]/60, mean_host, err, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'LineStyle','None',...
    'Color','b','LineWidth',1); 

mean_host_nophage = mean(cba38_nophage,2); 
std_host_nophage = std(cba38_nophage,0,2);

errorbar([0 175 945]/60, mean_host_nophage, std_host_nophage, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'b', ...
    'LineStyle','none', ...
    'Color','b','LineWidth',1); 


box on; xticks([0 2 4 6 8 10 12 14 16]);
yticks([1e-6 1e-3 1e0 1e3 1e6 1e9]);

%%  5) HP1 – H100

qpcr_host = 1e3*[1.03E+04	1.42E+04	1.78E+04;
8.99E+04	8.52E+04	1.30E+05;
9.45E+03	6.21E+03	6.66E+03];


qpcr_virus = 1e3*[3.03E+02	1.94E+02	4.94E+02;
1.38E+04	7.82E+03	1.07E+04;
8.77E+07	7.95E+07	8.48E+07];

moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;


load('./../data_2024/HP1_H100_2024.mat');

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
        host_den_seivd_paiwise_temp(host_den_seivd_paiwise_temp < cutoff_den) = 0;
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end


subplot(4,4,10)

for i = 1:length(Dc_array)
    hold on;
    if i == length(Dc_array)
        plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
    else
        plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:));
    end
end

hold on;
set(gca,'YScale','log');
xlabel('Time (hrs)')
ylabel('Phage/ml')
%title('PSA-HP1 – PSA H100');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
xlim([0 16]);

mean_virus = mean(qpcr_virus,2);             % mean across replicates (rows are time points)
err = std(qpcr_virus,0,2);          % std across replicates
errorbar([0 175 945]/60, mean_virus, err, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'LineStyle','none', ...
    'Color','r','LineWidth',1); 

ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]); xticks([0,2,4,6,8,10,12,14,16]);
box on; xticks([0 2 4 6 8 10 12 14 16])


subplot(4,4,9)
for i = 1:length(Dc_array)
    hold on;
    if i == length(Dc_array)
     

        plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
    else
        plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:));
    end
end



hold on;
set(gca,'YScale','log');
xlabel('Time (hrs)');
ylabel('Host/ml');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
%title('PSA-HP1 – PSA H100');
xlim([0 16]);

mean_host = mean(qpcr_host,2);             % mean across replicates (rows are time points)
err = std(qpcr_host,0,2);          % std across replicates
errorbar([0 175 945]/60, mean_host, err, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'LineStyle','None',...
    'Color','b','LineWidth',1); 

mean_host_nophage = mean(h100_nophage,2); 
std_host_nophage = std(h100_nophage,0,2);

errorbar([0 175 945]/60, mean_host_nophage, std_host_nophage, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'b', ...
    'LineStyle','none', ...
    'Color','b','LineWidth',1); 


box on; xticks([0 2 4 6 8 10 12 14 16]);yticks([1e-6 1e-3 1e0 1e3 1e6 1e9]);

%% 6) HP1 – 13-15


qpcr_host = 1e3*[1.48E+04	1.71E+04	2.01E+04
3.89E+04	4.62E+04	4.67E+04
7.19E+04	1.22E+05	8.05E+04];


qpcr_virus = 1e3*[4.44E+02	2.82E+02	4.34E+02;
1.07E+04	9.97E+03	1.61E+04;
5.53E+07	5.50E+07	4.84E+07];

moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;


load('./../data_2024/HP1_13-15_2024.mat');

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
        host_den_seivd_paiwise_temp(host_den_seivd_paiwise_temp < cutoff_den) = 0;
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end


subplot(4,4,12)
for i = 1:length(Dc_array)
    hold on;
    if i == length(Dc_array)
        plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
    else
        plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:));
    end
end


hold on;
set(gca,'YScale','log');
xlabel('Time (hrs)')
ylabel('Phage/ml')
%title('PSA-HP1 – PSA 13-15');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
xlim([0 16]);

mean_virus = mean(qpcr_virus,2);             % mean across replicates (rows are time points)
err = std(qpcr_virus,0,2);          % std across replicates
errorbar([0 175 945]/60, mean_virus, err, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'LineStyle','none', ...
    'Color','r','LineWidth',1); 

ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]); xticks([0,2,4,6,8,10,12,14,16]);
box on; xticks([0 2 4 6 8 10 12 14 16])


subplot(4,4,11)
for i = 1:length(Dc_array)
    hold on;
    if i == length(Dc_array)
     

        plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
    else
        plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:));
    end
end



hold on;
set(gca,'YScale','log');
xlabel('Time (hrs)');
ylabel('Host/ml');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
%title('PSA HP1 – PSA 13-15');
xlim([0 16]);
mean_host = mean(qpcr_host,2);             % mean across replicates (rows are time points)
err = std(qpcr_host,0,2);          % std across replicates
errorbar([0 175 945]/60, mean_host, err, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'LineStyle','None',...
    'Color','b','LineWidth',1); 

mean_host_nophage = mean(psa1315_nophage,2); 
std_host_nophage = std(psa1315_nophage,0,2);

errorbar([0 175 945]/60, mean_host_nophage, std_host_nophage, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'b', ...
    'LineStyle','none', ...
    'Color','b','LineWidth',1); 



box on; xticks([0 2 4 6 8 10 12 14 16]);yticks([1e-6 1e-3 1e0 1e3 1e6 1e9]);

%%  7) HS6 – H100

qpcr_host = 1e3*[1.33E+04	1.80E+04	1.69E+04;
7.74E+04	8.50E+04	8.66E+04;
2.35E+05	3.57E+05	2.43E+05];

qpcr_virus = 1e3*[7.79E+02	5.53E+02	7.96E+02;
3.95E+04	2.52E+04	1.86E+04;
1.12E+07	1.10E+07	9.59E+06];

moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;


load('./../data_2024/HS6_H100_2024.mat');

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
        host_den_seivd_paiwise_temp(host_den_seivd_paiwise_temp < cutoff_den) = 0;
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end


subplot(4,4,14)
for i = 1:length(Dc_array)
    hold on;
    if i == length(Dc_array)
        plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
    else
        plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:));
    end
end


hold on;
set(gca,'YScale','log');
xlabel('Time (hrs)')
ylabel('Phage/ml')
%title('PSA-HS6 – PSA H100');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
xlim([0 16]);

mean_virus = mean(qpcr_virus,2);             % mean across replicates (rows are time points)
err = std(qpcr_virus,0,2);          % std across replicates
errorbar([0 175 945]/60, mean_virus, err, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'LineStyle','none', ...
    'Color','r','LineWidth',1); 

ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]); xticks([0,2,4,6,8,10,12,14,16]);
box on; xticks([0 2 4 6 8 10 12 14 16])


subplot(4,4,13)
for i = 1:length(Dc_array)
    hold on;
    if i == length(Dc_array)
     

        plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
    else
        plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
            'LineWidth', 4, 'Color', gray(i,:));
    end
end


hold on;
set(gca,'YScale','log');
xlabel('Time (hrs)');
ylabel('Host/ml');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
%title('HS6 – H100');
xlim([0 16]);
mean_host = mean(qpcr_host,2);             % mean across replicates (rows are time points)
err = std(qpcr_host,0,2);          % std across replicates
errorbar([0 175 945]/60, mean_host, err, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k', ...
    'LineStyle','None',...
    'Color','b','LineWidth',1); 

mean_host_nophage = mean(h100_nophage,2); 
std_host_nophage = std(h100_nophage,0,2);

errorbar([0 175 945]/60, mean_host_nophage, std_host_nophage, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'b', ...
    'LineStyle','none', ...
    'Color','b','LineWidth',1); 


box on; xticks([0 2 4 6 8 10 12 14 16]);yticks([1e-6 1e-3 1e0 1e3 1e6 1e9]);

%%  8


qpcr_host = 1e3*[1.65E+04	2.02E+04	1.99E+04;
3.29E+04	3.08E+04	3.33E+04;
1.76E+05	1.48E+05	1.37E+05];

qpcr_virus = 1e3*[7.06E+02	7.25E+02	6.49E+02;
2.58E+04	2.93E+04	2.65E+04;
1.07E+07	1.25E+07	1.32E+07];


moi = mean(qpcr_virus(1,:))/mean(qpcr_host(1,:));
moi_list(end+1) = moi;


load('./../data_2024/HS6_13-15_2024.mat');


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
        host_den_seivd_paiwise_temp(host_den_seivd_paiwise_temp < cutoff_den) = 0;
    total_virus_seivd_pairwise{i} = total_virus_seivd_pairwise_temp;
    host_den_seivd_paiwise{i} = host_den_seivd_paiwise_temp;
    time_seivd{i} = time_seivd_temp;
end


% --- subplot 16: virus + its errorbars ---
subplot(4,4,16)
hold on
hLines16 = gobjects(length(Dc_array),1);
for i = 1:length(Dc_array)
    ls = '-'; if i==length(Dc_array), ls='--'; end
    hLines16(i) = plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
        'LineWidth',4,'Color',gray(i,:),'LineStyle',ls);
end
set(gca,'YScale','log'); xlim([0 16])
xlabel('Time (hrs)'); ylabel('Phage/ml'); %title('PSA-HS6 – PSA 13-15')
set(gca,'FontSize',20,'FontName','Times New Roman')
ylim([1e4 1e11]); yticks([1e4 1e6 1e8 1e10]); xticks(0:2:16); box on

mean_virus = mean(qpcr_virus,2);
err_v = std(qpcr_virus,0,2);
hPhage = errorbar([0 175 945]/60, mean_virus, err_v, ...
    's', 'MarkerSize',10, 'MarkerEdgeColor','k','MarkerFaceColor','k', ...
    'LineStyle','none', 'Color','r', 'LineWidth',1);   % red error bars

% --- subplot 15: host + its errorbars and legend ---
subplot(4,4,15)
hold on
hLines15 = gobjects(length(Dc_array),1);
for i = 1:length(Dc_array)
    ls = '-'; if i==length(Dc_array), ls='--'; end
    hLines15(i) = plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
        'LineWidth',4,'Color',gray(i,:),'LineStyle',ls);
end
set(gca,'YScale','log'); xlim([0 16])
xlabel('Time (hrs)'); ylabel('Host/ml'); %title('PSA-HS6 – PSA 13-15')
set(gca,'FontSize',20,'FontName','Times New Roman')
box on; xticks(0:2:16); yticks([1e-6 1e-3 1e0 1e3 1e6 1e9])

mean_host = mean(qpcr_host,2);
err_h = std(qpcr_host,0,2);
hHost = errorbar([0 175 945]/60, mean_host, err_h, ...
    's', 'MarkerSize',10, 'MarkerEdgeColor','k','MarkerFaceColor','k', ...
    'LineStyle','none', 'Color','b', 'LineWidth',1);   % blue error bars

mean_host_nophage = mean(h100_nophage,2); 
std_host_nophage = std(h100_nophage,0,2);

unHost = errorbar([0 175 945]/60, mean_host_nophage, std_host_nophage, ...
    's-', ...                      % square marker + line
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'b', ...
    'LineStyle','none', ...
    'Color','b','LineWidth',1); 


% Legend: lines + markers
% Dummy handles for headings
hModel = plot(nan,nan,'LineStyle','none','Marker','none'); % heading handle
hData = plot(nan,nan,'LineStyle','none','Marker','none'); % heading handle

% Collect handles
legHandles = [hModel; hLines15(:); hData; hPhage; hHost; unHost];

% Labels with headings
legLabels = { ...
    '\bf Model simulation', ...   % heading
    'Pairwise SEIVD; D_c = 10^5 cells/ml', ...
    'Pairwise SEIVD; D_c = 10^6 cells/ml', ...
    'Pairwise SEIVD; D_c = 5\times10^6 cells/ml', ...
    'Pairwise SEIV; D_c \rightarrow \infty', ...
    '\bf Data', ...               % heading
    'Phage density (qPCR-EXT)', ...
    'Infected Host density (qPCR-INT)',...
    'Uninfected Host density (qPCR-INT)'};

legend(legHandles, legLabels, 'Location','best', 'Box','on');


% subplot(4,4,16)
% for i = 1:length(Dc_array)
%     hold on;
%     if i == length(Dc_array)
%         plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
%             'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
%     else
%         plot(time_seivd{i}, total_virus_seivd_pairwise{i}, ...
%             'LineWidth', 4, 'Color', gray(i,:));
%     end
% end
% 
% 
% hold on;
% set(gca,'YScale','log');
% xlabel('Time (hrs)')
% ylabel('Phage/ml')
% %title('PSA-HS6 – PSA 13-15');
% set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
% xlim([0 16]);
% 
% mean_virus = mean(qpcr_virus,2);             % mean across replicates (rows are time points)
% err = std(qpcr_virus,0,2);          % std across replicates
% errorbar([0 175 945]/60, mean_virus, err, ...
%     's-', ...                      % square marker + line
%     'MarkerSize', 10, ...
%     'MarkerEdgeColor', 'k', ...
%     'MarkerFaceColor', 'k', ...
%     'LineStyle','none', ...
%     'Color','r','LineWidth',1); 
% 
% ylim([1e4 1e11]);
% yticks([1e4, 1e6, 1e8, 1e10]); xticks([0,2,4,6,8,10,12,14,16]);
% box on; xticks([0 2 4 6 8 10 12 14 16])
% 
% subplot(4,4,15)
% for i = 1:length(Dc_array)
%     hold on;
%     if i == length(Dc_array)
% 
% 
%         plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
%             'LineWidth', 4, 'Color', gray(i,:), 'LineStyle','--');
%     else
%         plot(time_seivd{i}, host_den_seivd_paiwise{i}, ...
%             'LineWidth', 4, 'Color', gray(i,:));
%     end
% end
% 
% 
% hold on;
% set(gca,'YScale','log');
% xlabel('Time (hrs)');
% ylabel('Host/ml');
% set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
% %title('PSA-HS6 – PSA 13-15');
% xlim([0 16]);
% mean_host = mean(qpcr_host,2);             % mean across replicates (rows are time points)
% err = std(qpcr_host,0,2);          % std across replicates
% errorbar([0 175 945]/60, mean_host, err, ...
%     's-', ...                      % square marker + line
%     'MarkerSize', 10, ...
%     'MarkerEdgeColor', 'k', ...
%     'MarkerFaceColor', 'k', ...
%     'LineStyle','None',...
%     'Color','b','LineWidth',1); 
% box on; xticks([0 2 4 6 8 10 12 14 16]);yticks([1e-6 1e-3 1e0 1e3 1e6 1e9]);
% 
% 
% % Display the legend for this subplot
% legend('D_c = 10^5 lysed cells/ml', ...
%      'D_c = 10^6 lysed cells/ml', ...
%       'D_c = 5\times10^6 lysed cells/ml', ...
%         ' SEIV (D_c \rightarrow \infty)')
%       % 'D_c = 10^7 lysed cells/ml', ...
% 
% 
% 
% legend('show');

%%

moi_list

moi_list.*infectious_ratio


%% 

drawnow;
set(findall(fig,'-property','FaceAlpha'),'FaceAlpha',1);
set(findall(fig,'-property','EdgeAlpha'),'EdgeAlpha',1);
set(findall(fig,'-property','AlphaData'),'AlphaData',[]);
