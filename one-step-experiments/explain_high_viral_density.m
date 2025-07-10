
clc;
clear all;
addpath(genpath(pwd)); % add current directory to path (mcmcstat is included)
addpath('./simulator/');
addpath('./mcmcstat/');

%%

B0 = 3e5;
V0 = 3e6*0.1;
theta_optimized = [ 0.2     2.18e-08         1.85        194.21          135];
infectious_ratio = 1;
NE_optimal = theta_optimized(5);
Dc = 5e11; %off
tvec = 0:0.1:15.8;

[time_seivd_temp, y_series_seivd] = one_step_simulate_noninfectious_seivd(tvec,B0,V0,theta_optimized,NE_optimal,Dc,infectious_ratio);
total_virus_seivd_pairwise_temp = sum(y_series_seivd(NE_optimal+3:NE_optimal+4,:));
host_den_seivd_paiwise_temp = sum(y_series_seivd(1:NE_optimal+2,:));
  
subplot(2,1,1)
plot(time_seivd_temp,host_den_seivd_paiwise_temp);
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);

  
subplot(2,1,2)
plot(time_seivd_temp,total_virus_seivd_pairwise_temp);
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage density (/ml)');
set(gca,'FontSize',20);
