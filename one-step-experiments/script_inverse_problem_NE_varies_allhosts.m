clc;
clear all;
addpath(genpath(pwd)); % add current directory to path (mcmcstat is included)
addpath('./simulator/');
addpath('/Users/rdey33/Downloads/MATLAB_DRIVE/mcmcstat/mcmcstat');

load("./results_same_phi/HP1_allhosts.mat");

load('parameters.mat','pars');


%% collecting statistics

free_phages_mean = mean(free_phages,2);
seed = 10;


data.ydata = free_phages_mean;
data.xdata = time_free_phages;

%% functions inclusion
%NE =85; %number of exposed class

moi_mean = mean(moi);
S0 = 1e8;
V0 = S0*moi_mean;



% theta given by OSU lab
beta_osu = 300;
r_osu = 0.265;
phi_osu = 1e-7;
tau_osu = 150  ;

theta_osu = [r_osu,phi_osu,tau_osu,beta_osu];


%% inference mcmc

model.ssfun = @(theta,data)  error_function_NE_varies(theta,data,S0,V0);


params = {
% initial values for the model states
    {'r', 0.24, 0, 0.5, 0.24,0.05}
    {'phi', 9.741e-8,  1e-10, 1e-6, 9.741e-8, 1e-7 }
    {'tau', 2.21,   0.25, 5, 2.21, 1 }
    {'beta', 61, 0, 700, 61, 100}
    {'NE',85,5,200,85,50};
    };




options.nsimu = 10000;
[results, chain] = mcmcrun(model,data,params,options);

chainstats(chain,results);



%%
format short E;
theta_inferred = median(chain(2000:end,:));

NE = round(theta_inferred(5));
y0(1) = S0;
y0(2:NE+2) = 0;
y0(NE+3) = V0;



dilution_factor = 100;
[time,y_series_inferred] = one_step_simulate(time_free_phages,y0,theta_inferred(1:4),NE,dilution_factor);
[time2,y_series_osu] = one_step_simulate(time_free_phages,y0,theta_osu,NE,dilution_factor);

%%
chain_effective = chain(2000:end,:);
for i=1:skips:length(chain_effective)

clear y
clear y0
NE =  round(chain_effective(i,5));
y(1,1) = S0;
y(1,2:NE+2) = 0;
y(1,NE+3) = V0;

[time2,y_series2,time_abs,pre_dil] = one_step_simulate(time_free_phages,y,chain_effective(i,:),NE, dilution_factor);

errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)
patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'linewidth',2,'edgealpha',0.2);hold on;


end

plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('PSA HP1--all hosts');
ylim([1e2 1e9]);
xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;