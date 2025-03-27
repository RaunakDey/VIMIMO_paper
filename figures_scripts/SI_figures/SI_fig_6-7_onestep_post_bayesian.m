%% Please ignore this script for now -- I will update this one later.

clc;
clear all;

%% trace, convergence tests, priors and posteriors.

%load('./results_same_phi/HS6_13-15_3-inferred.mat');
%13
%load("./../../one-step-experiments/result_replicates/18:2_18_12_100.mat");
%big problem
%load("./../../one-step-experiments/result_replicates/18:3_4_12.mat")
%load("./../../one-step-experiments/result_replicates/CBA18-3_18_2024_12.mat")
%load("./../../one-step-experiments/result_replicates/CBA38-1_38_13.mat")

%load("./../../one-step-experiments/result_replicates/HP1_H10011.mat");
%load("./../../one-step-experiments/result_replicates/HP1_13-1511.mat")
%load("./../../one-step-experiments/result_replicates/hs6_h100_18.mat")
%load("./../../one-step-experiments/result_replicates/HS6_13-1511.mat")



%% CBA 18:2 on CBA 18

load("./../../one-step-experiments/result_replicates/18:2_18_12_100.mat");
Si_onestep_utils
sgtitle('CBA 18:2 -- CBA 18 ');
set(gca, 'FontSize', 14, 'FontName','Times')



%% CBA 18:3 on CBA 4 -- redo this one.

load("./../../one-step-experiments/result_replicates/18:3_4_12.mat");
Si_onestep_utils
sgtitle('CBA 18:3 -- CBA 4 ');
set(gca, 'FontSize', 14, 'FontName','Times')



%% CBA 18:3 on CBA 18

load("./../../one-step-experiments/result_replicates/CBA38-1_38_13.mat");
Si_onestep_utils
sgtitle('CBA 38:1 -- CBA 38');
set(gca, 'FontSize', 14, 'FontName','Times')



%% CBA 38:1 on CBA 38


load("./../../one-step-experiments/result_replicates/CBA18-3_18_2024_12.mat");
Si_onestep_utils
sgtitle('CBA 18:3 -- CBA 18');
set(gca, 'FontSize', 14, 'FontName','Times')



%% HP1 on H100

load("./../../one-step-experiments/result_replicates/HP1_H10011.mat");
Si_onestep_utils
sgtitle('PSA HP1 -- PSA H100');
set(gca, 'FontSize', 14, 'FontName','Times')


%% HP1 on 13-15

load("./../../one-step-experiments/result_replicates/HP1_13-1511.mat");
Si_onestep_utils;
sgtitle('PSA HP1 -- PSA 13-15');
set(gca, 'FontSize', 14, 'FontName','Times')



%% HS6 on H100

load("./../../one-step-experiments/result_replicates/hs6_h100_18.mat");
Si_onestep_utils;
sgtitle('PSA HS6 -- PSA H100');
set(gca, 'FontSize', 14, 'FontName','Times')


%% HS6 on 13-15

load("./../../one-step-experiments/result_replicates/HS6_13-1511.mat");
Si_onestep_utils;
sgtitle('PSA HS6 -- PSA 13-15');
set(gca, 'FontSize', 14, 'FontName','Times')


