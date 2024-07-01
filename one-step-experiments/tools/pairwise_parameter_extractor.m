clear all;
clc;

r_bayesian = [];
phi_bayesian = [];
tau_bayesian = [];
beta_bayesian = [];
NE_bayesian = [];

r_bayesian_std = [];
phi_bayesian_std = [];
tau_bayesian_std = [];
beta_bayesian_std = [];
NE_bayesian_std = [];


load("./result_replicates/18:2_18_12_13.mat");
chain_effective = chain(5000:end,:);

r_bayesian(end+1) = median(chain_effective(:,1));
phi_bayesian(end+1) = median(chain_effective(:,2));
tau_bayesian(end+1) = median(chain_effective(:,3));
beta_bayesian(end+1) = median(chain_effective(:,4));
NE_bayesian(end+1) = median(chain_effective(:,5));

r_bayesian_std(end+1) = std(chain_effective(:,1));
phi_bayesian_std(end+1) = std(chain_effective(:,2));
tau_bayesian_std(end+1) = std(chain_effective(:,3));
beta_bayesian_std(end+1) = std(chain_effective(:,4));
NE_bayesian_std(end+1) = std(chain_effective(:,5));


load("./result_replicates/18:3_4_12.mat")
chain_effective = chain(5000:end,:);

r_bayesian(end+1) = median(chain_effective(:,1));
phi_bayesian(end+1) = median(chain_effective(:,2));
tau_bayesian(end+1) = median(chain_effective(:,3));
beta_bayesian(end+1) = median(chain_effective(:,4));
NE_bayesian(end+1) = median(chain_effective(:,5));

r_bayesian_std(end+1) = std(chain_effective(:,1));
phi_bayesian_std(end+1) = std(chain_effective(:,2));
tau_bayesian_std(end+1) = std(chain_effective(:,3));
beta_bayesian_std(end+1) = std(chain_effective(:,4));
NE_bayesian_std(end+1) = std(chain_effective(:,5));



load("./result_replicates/CBA18-3_18_2024_12.mat")
chain_effective = chain(5000:end,:);

r_bayesian(end+1) = median(chain_effective(:,1));
phi_bayesian(end+1) = median(chain_effective(:,2));
tau_bayesian(end+1) = median(chain_effective(:,3));
beta_bayesian(end+1) = median(chain_effective(:,4));
NE_bayesian(end+1) = median(chain_effective(:,5));

r_bayesian_std(end+1) = std(chain_effective(:,1));
phi_bayesian_std(end+1) = std(chain_effective(:,2));
tau_bayesian_std(end+1) = std(chain_effective(:,3));
beta_bayesian_std(end+1) = std(chain_effective(:,4));
NE_bayesian_std(end+1) = std(chain_effective(:,5));





load("./result_replicates/CBA38-1_38_13.mat")
chain_effective = chain(5000:end,:);

r_bayesian(end+1) = median(chain_effective(:,1));
phi_bayesian(end+1) = median(chain_effective(:,2));
tau_bayesian(end+1) = median(chain_effective(:,3));
beta_bayesian(end+1) = median(chain_effective(:,4));
NE_bayesian(end+1) = median(chain_effective(:,5));

r_bayesian_std(end+1) = std(chain_effective(:,1));
phi_bayesian_std(end+1) = std(chain_effective(:,2));
tau_bayesian_std(end+1) = std(chain_effective(:,3));
beta_bayesian_std(end+1) = std(chain_effective(:,4));
NE_bayesian_std(end+1) = std(chain_effective(:,5));




load("./result_replicates/HP1_H10011.mat");
chain_effective = chain(5000:end,:);

r_bayesian(end+1) = median(chain_effective(:,1));
phi_bayesian(end+1) = median(chain_effective(:,2));
tau_bayesian(end+1) = median(chain_effective(:,3));
beta_bayesian(end+1) = median(chain_effective(:,4));
NE_bayesian(end+1) = median(chain_effective(:,5));

r_bayesian_std(end+1) = std(chain_effective(:,1));
phi_bayesian_std(end+1) = std(chain_effective(:,2));
tau_bayesian_std(end+1) = std(chain_effective(:,3));
beta_bayesian_std(end+1) = std(chain_effective(:,4));
NE_bayesian_std(end+1) = std(chain_effective(:,5));




load("./result_replicates/HP1_13-1511.mat")
chain_effective = chain(5000:end,:);

r_bayesian(end+1) = median(chain_effective(:,1));
phi_bayesian(end+1) = median(chain_effective(:,2));
tau_bayesian(end+1) = median(chain_effective(:,3));
beta_bayesian(end+1) = median(chain_effective(:,4));
NE_bayesian(end+1) = median(chain_effective(:,5));

r_bayesian_std(end+1) = std(chain_effective(:,1));
phi_bayesian_std(end+1) = std(chain_effective(:,2));
tau_bayesian_std(end+1) = std(chain_effective(:,3));
beta_bayesian_std(end+1) = std(chain_effective(:,4));
NE_bayesian_std(end+1) = std(chain_effective(:,5));



%load("./result_replicates/HS6_H100_11.mat")
load("./result_replicates/hs6_h100_15.mat");
chain_effective = chain(5000:end,:);

r_bayesian(end+1) = median(chain_effective(:,1));
phi_bayesian(end+1) = median(chain_effective(:,2));
tau_bayesian(end+1) = median(chain_effective(:,3));
beta_bayesian(end+1) = median(chain_effective(:,4));
NE_bayesian(end+1) = median(chain_effective(:,5));

r_bayesian_std(end+1) = std(chain_effective(:,1));
phi_bayesian_std(end+1) = std(chain_effective(:,2));
tau_bayesian_std(end+1) = std(chain_effective(:,3));
beta_bayesian_std(end+1) = std(chain_effective(:,4));
NE_bayesian_std(end+1) = std(chain_effective(:,5));



load("./result_replicates/HS6_13-1511.mat")
chain_effective = chain(5000:end,:);

r_bayesian(end+1) = median(chain_effective(:,1));
phi_bayesian(end+1) = median(chain_effective(:,2));
tau_bayesian(end+1) = median(chain_effective(:,3));
beta_bayesian(end+1) = median(chain_effective(:,4));
NE_bayesian(end+1) = median(chain_effective(:,5));

r_bayesian_std(end+1) = std(chain_effective(:,1));
phi_bayesian_std(end+1) = std(chain_effective(:,2));
tau_bayesian_std(end+1) = std(chain_effective(:,3));
beta_bayesian_std(end+1) = std(chain_effective(:,4));
NE_bayesian_std(end+1) = std(chain_effective(:,5));

