clc;
clear all;

addpath(genpath('./mcmcstat/'));
addpath('./revision/corrected4/');

%Bayesian data file.
load('./revision/corrected4/v23.mat');
%%

%[mcmcresults, chain, s2chain]= mcmcrun(mcmcmodel,data,mcmcparam,mcmcoptions);

beta_lower = 1;
phi_lower = -9;
tau_lower = 0.5;
r_lower = 0.1;
dc_lower = 6;
parameter_lower_bound = [repmat(beta_lower ,1,9), repmat(phi_lower ,1,9), repmat(tau_lower ,1,9), repmat(r_lower ,1,5),repmat(dc_lower ,1,5) ];


beta_upper = 700;
phi_upper = -5;
tau_upper = 8;
r_upper = 1;
dc_upper = 8;
parameter_upper_bound = [repmat(beta_upper ,1,9), repmat(phi_upper ,1,9), repmat(tau_upper ,1,9), repmat(r_upper ,1,5),repmat(dc_upper ,1,5)];

beta_std = 500;
phi_std = 2;
tau_std = 5;
r_std = 0.05;
dc_std = 1;

theta_std = [repmat(beta_std ,1,9), repmat(phi_std ,1,9), repmat(tau_std ,1,9), repmat(r_std ,1,5),repmat(dc_std ,1,5)];


%standard normal lower and upper
sn_lower = (parameter_lower_bound  - theta_start)./theta_std;
sn_upper = (parameter_upper_bound  - theta_start)./theta_std;

% z = (x - mu)./sigma

for i = 1:37
    mcmcparam{1,i}{1,3} = sn_lower(i);
    mcmcparam{1,i}{1,4} = sn_upper(i);
end

%%

mcmcmodel.sigma2 = 150;
mcmcmodel.S20 = 150;
mcmcmodel.N0 = 1; % noninformative

mcmcoptions.nsimu = 10000;
[mcmcresults_continued, chain_continues, s2chain_continues]= mcmcrun(mcmcmodel,data,mcmcparam,mcmcoptions);


%% sanity checks -- only for Raunak -- do not run

% ## How I did the Bayesian inference
% The parameters from One step has so much deviation from the likelihood
% curve that it will lead to a flat posterior, extremely low acceptance
% rates.

% Used gradient descent on one of the datasets to fit parameters. 



%% Model and inference settings

% 1. Differential equations: written in the paper, used model =
% SEIVD_diff_NE_diff_debris_abs 
% used NE = 200 for all cases.
% Did use same \beta for the whole of the differential equations.
% 2. Bayesian model: 
% set sigma, S0, N0 as hyperpriors of the priors of the sigmaLL (equal for
% all datapoints).
% 3. Inference settings: dram method, set adaptive steps to 100.
% these settings are already saved.
