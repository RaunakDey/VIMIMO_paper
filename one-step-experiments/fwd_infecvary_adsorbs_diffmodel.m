clc;
clear all;
addpath(genpath(pwd)); % add current directory to path (mcmcstat is included)
addpath('./simulator/');
addpath('./mcmcstat/');


%%%% just looking at one phage for now please.

load('./data_2024/CBA18-2_18_2024.mat');
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');
num_replicates = 1;


% parameters.
r_optimal = 0.245;
phi_optimal = 2.18e-8;
tau_optimal = 1.85;
beta_optimal = 94.2;
NE_optimal = 135;
prob_effective_infection = 0.5;
prob_lysis_reset = 0e-7;



theta_optimized(6) = prob_effective_infection; %prob effective inhibition
theta_optimized(7) = prob_lysis_reset; %prob lysis reset

NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);


dilution_factor = 100;


tvec = 0:0.1:50;
[time_opt,y_series_opt] = one_step_simulate_infec_vary_adsorbs_diffmodel(tvec,y0,theta_optimized,NE_optimal,dilution_factor);


total_virus = y_series_opt(end,:)+y_series_opt(end-1,:);

subplot(2,1,1)

plot(time_opt,total_virus,'-g','LineWidth',2);hold on;
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
xlim([0 tvec(end)]);

ylim([1e3 1e9]);
yticks([1e3 1e5 1e7 1e9])


host_den = sum(y_series_opt(1:end-2,:));

subplot(2,1,2)
plot(time_opt,host_den,'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
xlabel('Time (hr)');
ylabel('Host density (/ml)');
set(gca,'FontSize',20);
title('CBA18-2 on CBA 18');
xlim([0 tvec(end)]);
ylim([1e-1 1e7]);
yticks([1e-1, 1e1, 1e3, 1e5, 1e7]);

