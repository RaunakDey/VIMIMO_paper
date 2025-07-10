clc;
clear all;
addpath(genpath(pwd)); % add current directory to path (mcmcstat is included)
addpath('./simulator/');
addpath('./mcmcstat/');
cutoff = 1e-2;

%%
MOI_list = [10,1,0.1];
color_list = ['r','b','k'];

for iter = 1:length(MOI_list)

MOI  = MOI_list(iter);
color_line = color_list(iter);

B0 = 2e6;
V0 = B0*MOI;
theta_optimized = [ 0.2     2.18e-08         1.85        194.21          135];
NE_optimal = theta_optimized(5);

tvec = 0:0.1:15.8;

[time_seiv, y_series_seiv] = simulate_pairwise_seiv(tvec,B0,V0,theta_optimized,NE_optimal);

virus_den_seiv = y_series_seiv(NE_optimal+3,:);
host_den_seiv = sum(y_series_seiv(1:NE_optimal+2,:));
host_den_seiv(host_den_seiv<cutoff) = 0;
  

subplot(1,2,1)
plot(time_seiv,host_den_seiv, Color=color_line)
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
hold on;

  
subplot(1,2,2)
plot(time_seiv,virus_den_seiv, Color=color_line);
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Phage density (/ml)');
set(gca,'FontSize',20);
hold on;

end



