clc;
clear all;

addpath(genpath('./mcmcstat/'));
addpath('./revision/corrected4/');

%Bayesian data file.
load('./revision/corrected4/v23.mat');


%[mcmcresults, chain, s2chain]= mcmcrun(mcmcmodel,data,mcmcparam,mcmcoptions);
mcmcoptions.nsimu = 1000;
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
