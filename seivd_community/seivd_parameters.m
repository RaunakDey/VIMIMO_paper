clear all;
close all;
clc;


%%
load('./../community/data/qpcr.mat');
load('./../community/results/v23.mat');

model =  SEIVD_diff_NE_diff_debris_abs(5,5,200);
model.name = 'SEIVD-diffabs';
model.debris_inhib = 2;
model.debris_inhib2 = 2;
model.debris_inhib3 = 2;
model.debris_inhib4 = 2;
model.debris_inhib5 = 2;



tvec = 0:0.05:15.75; % for better viz

sd_theta = [10.0000
    0.5000
   10.0000
   10.0000
    5.0000
   10.0000
   10.0000
   10.0000
   10.0000
    0.2000
    0.2000
    0.2000
    0.2000
    0.2000
    0.2000
    0.2000
    0.2000
    0.2000
    0.1000
    0.1000
    0.1000
    0.1000
    0.1000
    0.1000
    0.1000
    0.1000
    0.1000
    0.0500
    0.0500
    0.0500
    0.0500
    0.0500
    0.1000
    0.1000
    0.1000
    0.1000
    0.1000];

num_steps = 30;
step_size = 1;




theta_optimized = [220.3648
    0.3755
  166.1027
   79.6696
   19.9832
  481.9219
  534.6305
   30.6009
   51.4209
   -7.7640
   -7.2083
   -7.1083
   -7.6198
   -7.0464
   -7.2306
   -7.2277
   -7.8731
   -7.6460
    2.6021
    3.0896
    3.1028
    2.9918
    1.8965
    1.7983
    3.2942
    5.2070
    2.0063
    0.1760
    0.1706
    0.2450
    0.6135
    0.5282
    6.8018
    6.8782
    7.0669
    6.2867
    6.2053];

theta_initial = theta_optimized;

theta_optimized = search_minimum_community(theta_initial,sd_theta,data,model,pars2,mcmcpars,num_steps, step_size);





%% Settings for running the model



% the important parameters.
pars_afterinf.NH = 5;
pars_afterinf.NV = 5;
pars_afterinf.r = [0.17692
      0.22069
      0.29393
      0.66577
      0.52807];

pars_afterinf.M =[0   1   0   0   0
   1   1   1   0   0
   0   0   1   0   0
   0   0   0   1   1
   0   0   0   1   1];
pars_afterinf.NE = 200*pars_afterinf.M;

pars_afterinf.beta= [         0       2.8875            0            0            0
        194.9       204.58       100.46            0            0
            0            0       19.922            0            0
            0            0            0       525.38       60.653
            0            0            0       488.06       51.326];

pars_afterinf.tau = [  0       2.9978            0            0            0
       1.7394        2.746       2.3189            0            0
            0            0        1.988            0            0
            0            0            0        1.822       4.7139
            0            0            0       2.3157       1.9868];

% BIG PROBLEM HERE.
pars_afterinf.eta= [         0      0.33333            0            0            0
      0.58824      0.37037      0.43478            0            0
            0            0          0.5            0            0
            0            0            0      0.55556      0.21277
            0            0            0      0.43478          0.5];


pars_afterinf.phi = [            0   5.8924e-08            0            0            0
   1.5318e-08   7.7883e-08   2.4067e-08            0            0
            0            0   7.9035e-08            0            0
            0            0            0    6.114e-08   1.1905e-08
            0            0            0   6.0123e-08   2.2635e-08];

%pars_afterinf.phi = zeros(5,5);

pars_afterinf.Dc = 5.0415e+06;
pars_afterinf.Dc2 = 5.9627e+06;
pars_afterinf.Dc3 =  1.1927e+07;
pars_afterinf.Dc4 = 1.8324e+06;
pars_afterinf.Dc5 =  1.52e+06;


pars_afterinf.V0 = [4.2887e+05
   2.8689e+05
     5.28e+05
   1.1033e+05
    1.151e+07];

pars_afterinf.S0 = [2.5111e+06
   5.6423e+06
   3.0257e+06
    6.205e+06
   7.7533e+06];


pars_afterinf.V0 =  V0;
pars_afterinf.S0 =  S0;

pars_afterinf.epsilon = [1 1 1 1 1 1 1 1 1 1];

% not important -- the model class does not call these paramters, the
% plots should not change if these parameters are changed.
pars_afterinf.a = eye(5);
pars_afterinf.m = [  0.00081472
   0.00090579
   0.00012699
   0.00091338
   0.00063236];

pars_afterinf.m = zeros(5,1);
pars_afterinf.q = [0.5
          0.5
          0.5
          0.5
          0.5];

pars_afterinf.q = zeros(5,1);
pars_afterinf.prob = [    0
     0
     0
     0
     0];
