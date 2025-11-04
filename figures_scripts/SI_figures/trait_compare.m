clear all;
clc;

%% SEIVD
x=1:9;


load('./../../community/results/v25.mat');

% change the log transformed part of the chain to original
chain_comm = chain(30001:50000,1:27);
chain_comm(:,10:18) = 10.^chain_comm(:,10:18);



%% CBA 18 -- phi18:2


load("./../../one-step-experiments/result_replicates/18:2_18_12_100.mat");
chain_pairwise = chain(30001:end,:);


% Example with 9 side-by-side violin comparisons
figure(1)
t = tiledlayout(1,9,'TileSpacing','compact','Padding','compact');
ylabel(t,'Burst size');   % shared y label
set(gca,'FontSize',18,'FontName','Times')


figure(2)
t = tiledlayout(1,9,'TileSpacing','compact','Padding','compact');
ylabel(t,'\tau');   % shared y label
set(gca,'FontSize',18,'FontName','Times')

figure(3)
t = tiledlayout(1,9,'TileSpacing','compact','Padding','compact');
ylabel(t,'\phi');   % shared y label
set(gca,'FontSize',18,'FontName','Times')



% phi
chainA = chain_pairwise(:,2);
chainB = chain_comm(:,10);
fprintf('cba 18:2 on cba 18 -- phi\n');
res = bayes_eq_test_deterministic(chainA, chainB);
figure(3)
nexttile
%violin2chains(log(chainA)/2.303, log(chainB)/2.303);
violin2chains(chainA, chainB);
set(gca,'FontSize',18,'FontName','Times','YScale','log')
title('\phi18:2 – CBA 18');


% tau
chainA = chain_pairwise(:,3);
chainB = chain_comm(:,19);
fprintf('cba 18:2 on cba 18 -- tau\n');
res = bayes_eq_test_deterministic(chainA, chainB);

figure(2)
nexttile
violin2chains(chainA, chainB);
%ylim([0.5 5.5]);
set(gca,'FontSize',18,'FontName','Times')
title('\phi18:2 – CBA 18');


% beta

chainA = chain_pairwise(:,4);
chainB = chain_comm(:,1);
fprintf('cba 18:2 on cba 18 -- beta\n');
res = bayes_eq_test_deterministic(chainA, chainB);

figure(1)
nexttile
violin2chains(chainA, chainB);
%ylim([0 700]);
set(gca,'FontSize',18,'FontName','Times')
title('\phi18:2 – CBA 18');

%% CBA 4 -- phi18:3



%load("./../../one-step-experiments/result_replicates/18:3_4_104.mat");
load("./../../one-step-experiments/result_replicates/18:3_4_12.mat")
chain_pairwise = chain(10000:end,:);

% phi
chainA = chain_pairwise(:,2);
chainB = chain_comm(:,11);
fprintf('cba 18:3 on cba 4 -- phi\n');
res = bayes_eq_test_deterministic(chainA, chainB);

figure(3)
nexttile
violin2chains(chainA, chainB);
set(gca,'FontSize',18,'FontName','Times','YScale','log')

set(gca,'FontSize',18,'FontName','Times')
title('\phi18:2 – CBA 18');




% tau
chainA = chain_pairwise(:,3);
chainB = chain_comm(:,20);
fprintf('cba 18:3 on cba 4  -- tau\n');
res = bayes_eq_test_deterministic(chainA, chainB);
figure(2)
nexttile
violin2chains(chainA, chainB);
%ylim([0.5 5.5]);
set(gca,'FontSize',18,'FontName','Times')
title('\phi18:2 – CBA 18');


% beta
chainA = chain_pairwise(:,4);
chainB = chain_comm(:,2);
fprintf('cba 18:3 on cba 4  -- beta\n');
res = bayes_eq_test_deterministic(chainA, chainB);
figure(1)
nexttile
violin2chains(chainA, chainB);

set(gca,'FontSize',18,'FontName','Times')
title('\phi18:3 – CBA 4');
%% CBA 18 -- phi18:3

load("./../../one-step-experiments/result_replicates/CBA18-3_18_2024_12.mat")
chain_pairwise = chain(5001:end,:);


% phi
chainA = chain_pairwise(:,2);
chainB = chain_comm(:,12);
fprintf('cba 18:3 on cba 18 -- phi\n');
res = bayes_eq_test_deterministic(chainA, chainB);
figure(3)
nexttile
violin2chains(chainA, chainB);
set(gca,'FontSize',18,'FontName','Times','YScale','log')
title('\phi18:2 – CBA 18');



% tau
chainA = chain_pairwise(:,3);
chainB = chain_comm(:,21);
fprintf('cba 18:3 on cba 18  -- tau\n');
res = bayes_eq_test_deterministic(chainA, chainB);

figure(2)
nexttile
violin2chains(chainA, chainB);
ylim([0.5 5.5]);
set(gca,'FontSize',18,'FontName','Times')
title('\phi18:2 – CBA 18');


% beta
chainA = chain_pairwise(:,4);
chainB = chain_comm(:,3);
fprintf('cba 18:3 on cba 18  -- beta\n');
res = bayes_eq_test_deterministic(chainA, chainB);

figure(1)
nexttile
violin2chains(chainA, chainB);
ylim([0 700]);
set(gca,'FontSize',18,'FontName','Times')
title('\phi18:3 – CBA 18');
%% CBA 38 -- phi38:1

load("./../../one-step-experiments/result_replicates/CBA38-1_38_13.mat")
chain_pairwise = chain(5001:end,:);

% phi
chainA = chain_pairwise(:,2);
chainB = chain_comm(:,14);
fprintf('cba 38:1 on cba 38 -- phi\n');
res = bayes_eq_test_deterministic(chainA, chainB);
figure(3)
nexttile
violin2chains(chainA, chainB);
set(gca,'FontSize',18,'FontName','Times','YScale','log')
set(gca,'FontSize',18,'FontName','Times')
title('\phi38:1 – CBA 38');


% tau
chainA = chain_pairwise(:,3);
chainB = chain_comm(:,23);
fprintf('cba 38:1 on cba 38  -- tau\n');
res = bayes_eq_test_deterministic(chainA, chainB);
figure(2)
nexttile
violin2chains(chainA, chainB);
ylim([0.5 5.5]);
set(gca,'FontSize',18,'FontName','Times')
title('\phi38:1 – CBA 38');


% beta
chainA = chain_pairwise(:,4);
chainB = chain_comm(:,5);
fprintf('cba 38:1 on cba 38  -- beta\n');
res = bayes_eq_test_deterministic(chainA, chainB);

figure(1)
nexttile
violin2chains(chainA, chainB);
%ylim([0 700]);
set(gca,'FontSize',18,'FontName','Times')
title('\phi38:1 – CBA 38');
%% PSA HP1 on H100

load("./../../one-step-experiments/result_replicates/HP1_H10011.mat");
chain_pairwise = chain(5001:end,:);

% phi
chainA = chain_pairwise(:,2);
chainB = chain_comm(:,15);
fprintf('PSA HP1 on PSA H100 -- phi\n');
res = bayes_eq_test_deterministic(chainA, chainB);
figure(3)
nexttile
violin2chains(chainA, chainB);
set(gca,'FontSize',18,'FontName','Times','YScale','log')
title('PSA-HP1 – PSA H100');
set(gca,'YScale','log');


% tau
chainA = chain_pairwise(:,3);
chainB = chain_comm(:,24);
fprintf('PSA HP1 on PSA H100 -- tau\n');
res = bayes_eq_test_deterministic(chainA, chainB);
figure(2)
nexttile
violin2chains(chainA, chainB);
%ylim([0.5 5.5]);
set(gca,'FontSize',18,'FontName','Times')
title('PSA-HP1 – PSA H100');


% beta
chainA = chain_pairwise(:,4);
chainB = chain_comm(:,6);
fprintf('PSA HP1 on PSA H100  -- beta\n');
res = bayes_eq_test_deterministic(chainA, chainB);

figure(1)
nexttile
violin2chains(chainA, chainB);
%ylim([0 700]);
set(gca,'FontSize',18,'FontName','Times');
title('PSA-HP1 – PSA H100');
%% PSA HP1 on 13-15
load("./../../one-step-experiments/result_replicates/HP1_13-1511.mat")
chain_pairwise = chain(5001:end,:);

% phi
chainA = chain_pairwise(:,2);
chainB = chain_comm(:,16);
fprintf('PSA HP1 on PSA H100 -- phi\n');
res = bayes_eq_test_deterministic(chainA, chainB);
figure(3)
nexttile
violin2chains(chainA, chainB);
set(gca,'FontSize',18,'FontName','Times','YScale','log')
title('PSA-HP1 – PSA 13-15');
set(gca,'YScale','log');




% tau
chainA = chain_pairwise(:,3);
chainB = chain_comm(:,25);
fprintf('PSA HP1 on PSA H100 -- tau\n');
res = bayes_eq_test_deterministic(chainA, chainB);
figure(2)
nexttile
violin2chains(chainA, chainB);
%ylim([0.5 5.5]);
set(gca,'FontSize',18,'FontName','Times')
title('PSA-HP1 – PSA 13-15');


% beta
chainA = chain_pairwise(:,4);
chainB = chain_comm(:,7);
fprintf('PSA HP1 on PSA H100  -- beta\n');
res = bayes_eq_test_deterministic(chainA, chainB);

figure(1)
nexttile
violin2chains(chainA, chainB);
%ylim([0 700]);
set(gca,'FontSize',18,'FontName','Times')
title('PSA-HP1 – PSA 13-15');
%% PSA HS6 on H100


load("./../../one-step-experiments/result_replicates/hs6_h100_18.mat")
chain_pairwise = chain(5001:end,:);

% phi
chainA = chain_pairwise(:,2);
chainB = chain_comm(:,17);
fprintf('PSA HS6 on PSA H100 -- phi\n');
res = bayes_eq_test_deterministic(chainA, chainB);
figure(3)
nexttile
violin2chains(chainA, chainB);
set(gca,'FontSize',18,'FontName','Times','YScale','log')
title('PSA-HS6 – PSA H100');
set(gca,'YScale','log');




% tau
chainA = chain_pairwise(:,3);
chainB = chain_comm(:,26);
fprintf('PSA HS6 on PSA H100 -- tau\n');
res = bayes_eq_test_deterministic(chainA, chainB);
figure(2)
nexttile
violin2chains(chainA, chainB);
%ylim([0.5 5.5]);
set(gca,'FontSize',18,'FontName','Times')
title('PSA-HS6 – PSA H100');



% beta
chainA = chain_pairwise(:,4);
chainB = chain_comm(:,8);
fprintf('PSA HS6 on PSA H100  -- beta\n');
res = bayes_eq_test_deterministic(chainA, chainB);

figure(1)
nexttile
violin2chains(chainA, chainB);
%ylim([0 700]);
set(gca,'FontSize',18,'FontName','Times')
title('PSA-HS6 – PSA H100');
%% PSA HS6 on 13-15

load("./../../one-step-experiments/result_replicates/HS6_13-1511.mat");
chain_pairwise = chain(5001:end,:);

% phi
chainA = chain_pairwise(:,2);
chainB = chain_comm(:,18);
fprintf('PSA HS6 on PSA 13-15 -- phi\n');
res = bayes_eq_test_deterministic(chainA, chainB);
figure(3)
nexttile
violin2chains(chainA, chainB);
set(gca,'FontSize',18,'FontName','Times','YScale','log')
title('PSA-HS6 – PSA 13-15');
set(gca,'YScale','log');
set(gca,'YScale','log');




% tau
chainA = chain_pairwise(:,3);
chainB = chain_comm(:,27);
fprintf('PSA HS6 on PSA 13-15 -- tau\n');
res = bayes_eq_test_deterministic(chainA, chainB);
figure(2)
nexttile
violin2chains(chainA, chainB);
%ylim([0.5 5.5]);
set(gca,'FontSize',18,'FontName','Times')
title('PSA-HS6 – PSA 13-15');


% beta
chainA = chain_pairwise(:,4);
chainB = chain_comm(:,9);
fprintf('PSA HS6 on PSA 13-15  -- beta\n');
res = bayes_eq_test_deterministic(chainA, chainB);

figure(1)
nexttile
violin2chains(chainA, chainB);
%ylim([0 700]);
set(gca,'FontSize',18,'FontName','Times');
title('PSA-HS6 – PSA 13-15');

%% PSA additional HP1


load("./../../one-step-experiments/result_replicates/HP1_H10011.mat");
chain_pairwise = chain(5001:end,:);
chainA = chain_pairwise(:,4);
chainB = chain_comm(:,6);

load("./../../one-step-experiments/result_replicates/hp1_all13.mat");
hp1_all_beta = chain(5001:end,4);
median(hp1_all_beta)

res = bayes_eq_test_deterministic(chainA, hp1_all_beta,'mode','std')
res = bayes_eq_test_deterministic(chainB, hp1_all_beta,'mode','std')




load("./../../one-step-experiments/result_replicates/HP1_13-1511.mat")
chain_pairwise = chain(5001:end,:);
chainC = chain_pairwise(:,4);
chainD = chain_comm(:,7);


res = bayes_eq_test_deterministic(chainC, hp1_all_beta,'mode','std')
res = bayes_eq_test_deterministic(chainD, hp1_all_beta,'mode','std')



colors = [0.5 0.5 0.5;        % dark gray
          171/255 193/255 157/255; % greenish
          0.5 0.5 0.5;
          171/255 193/255 157/255; % greenish
          0.8 0.8 0.8;
                    0.8 0.8 0.8;];       % lighter gray
figure(5)
res = violinNchains({chainA, chainB,chainC, chainD, hp1_all_beta}, ...
              'labels',{'PSA H100 (in pairwise)','PSA H100 (in community)','PSA 13-15 (in pairwise)','PSA 13-15 (in community)','all hosts'}, ...
              'colors',colors, ...
              'width',0.25, 'alpha',0.8, ...
              'support',[0 Inf]);
set(gca,'FontSize',18,'FontName','Times')
xtickangle(90);
ylabel('Burst size (virions/cell)')

%% PSA additional HS6


load("./../../one-step-experiments/result_replicates/hs6_h100_18.mat")
chain_pairwise = chain(5001:end,:);
chainA = chain_pairwise(:,4);
chainB = chain_comm(:,9);


load("./../../one-step-experiments/result_replicates/hs6_all13.mat");
hs6_all_beta = chain(burn:end,4);
res = bayes_eq_test_deterministic(chainA, hs6_all_beta,'mode','std')
res = bayes_eq_test_deterministic(chainB, hs6_all_beta,'mode','std')
median(hs6_all_beta)

load("./../../one-step-experiments/result_replicates/HS6_13-1511.mat");
chain_pairwise = chain(5001:end,:);
chainC = chain_pairwise(:,4);
chainD = chain_comm(:,9);
res = bayes_eq_test_deterministic(chainC, hs6_all_beta, 'mode','std')
res = bayes_eq_test_deterministic(chainD, hs6_all_beta, 'mode','std')


colors = [0.5 0.5 0.5;        % dark gray
          171/255 193/255 157/255; % greenish
          0.5 0.5 0.5;
          171/255 193/255 157/255; % greenish
          0.8 0.8 0.8;
                    0.8 0.8 0.8;];       % lighter gray


figure(6)
res = violinNchains({chainA, chainB,chainC, chainD, hs6_all_beta}, ...
              'labels',{'PSA H100 (in pairwise)','PSA H100 (in community)','PSA 13-15 (in pairwise)','PSA 13-15 (in community)','all hosts'}, ...
              'colors',colors, ...
              'width',0.25, 'alpha',0.8, ...
              'support',[0 Inf]);
set(gca,'FontSize',18,'FontName','Times')
xtickangle(90);
ylabel('Burst size (virions/cell)')