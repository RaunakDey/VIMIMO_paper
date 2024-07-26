clc;
clear all;
addpath(genpath(pwd)); % add current directory to path (mcmcstat is included)
addpath('./simulator/');
addpath('./mcmcstat/');


%load('./data_2024/HP1_H100_2024.mat');
%load('./data_2024/HP1_13-15_2024.mat');
load('./data_2024/HS6_H100_2024.mat');
%load('./data_2024/HS6_13-15_2024.mat');


name = string(labels.phage)+'_'+string(labels.host);
time_free_phages = time_free_phages/60;
load('parameters.mat','pars');


format short g


seed = 18;

num_replicates = length(V0_replicates)/3;



%% plot optimal

NE_optimal = round(theta_optimized(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE_optimal+2) = 0;
y0(NE_optimal+3) = mean(V0_replicates);

dilution_factor = 100;

tvec = 0:0.01:3.5;
[time_opt,y_series_opt] = one_step_simulate(tvec,y0,theta_optimized,NE_optimal,dilution_factor);

figure(10)
plot(time_opt,y_series_opt(end,:),'-k','LineWidth',2);hold on;
set(gca,'YScale','log');
for i = 1:num_replicates
    for j = 1:3
    plot(data.xdata{i}./60,data.ydata{i}(:,j),'bo'); hold on;
    end
end

%% inference mcmc

model.ssfun = @(theta,data)  error_triplicates_multiples(theta,data,S0_replicates,V0_replicates);

r_opt = theta_optimized(1);
phi_opt = theta_optimized(2);
tau_opt = theta_optimized(3);
beta_opt = theta_optimized(4);
NE_opt = round(theta_optimized(5));

r_osu = theta_osu(1);
phi_osu = theta_osu(2);
tau_osu = theta_osu(3);
beta_osu = theta_osu(4);


params = {
% initial values for the model states
    {'r', r_opt, 0, 0.5, r_osu,0.05}
    {'phi', phi_opt,  1e-10, 1e-6,phi_osu, 1e-7 }
    {'tau', tau_opt,   0.25, 5, tau_osu, 5 }
    {'beta', beta_opt, 0, 700, 310, 150}
    {'NE',NE_opt,5,400,NE_opt,100};
    };


error_prior = error_triplicates_multiples(theta_optimized,data,S0_replicates,V0_replicates)

 
% model.S20 = 1;
% model.N0 = 4;
% options.updatesigma = 1;


model.S20 = 1;
model.S20 = error_prior;
model.N0 = 1;
options.updatesigma = 1;

%new
model.S20 = 1;
seed = 100;

options.nsimu =  50000;
[results, chain, s2chain] = mcmcrun(model,data,params,options);

chainstats(chain,results);


%% plots
burn = options.nsimu*0.5;
theta_inferred = median(chain(burn:end,:));

NE = round(theta_inferred(5));
clear y0
y0(1) = mean(S0_replicates);
y0(2:NE+2) = 0;
y0(NE+3) = mean(V0_replicates);




dilution_factor = 100;
[time,y_series_inferred] = one_step_simulate(time_free_phages,y0,theta_inferred(1:4),NE,dilution_factor);

%% saving

%save("./result_replicates/hs6_h100_18.mat");
save("./result_replicates/hs6_h100_"+string(seed));
