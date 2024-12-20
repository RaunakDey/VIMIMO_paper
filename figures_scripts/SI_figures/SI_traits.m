clear all;
clc;

%% SEIVD
x=1:9;



%load('./v14-4.mat');
load('./../../community/results/v23.mat');

% change the log transformed part of the chain to original
chain(:,10:18) = exp(2.303.*chain(:,10:18));

%start_new = 7000;
start_new = 30001;
beta(:,1) = median(chain(start_new:end,1:9)) ;
r(:,1) = median(chain(start_new:end,28:32)) ;
phi(:,1) = median(chain(start_new:end,10:18)) ;
tau(:,1) = median(chain(start_new:end,19:27)) ;

error_beta_seivd = std(chain(start_new:end,1:9)) ;
error_r_seivd = std(chain(start_new:end,28:32)) ;
error_phi_seivd = std(chain(start_new:end,10:18)) ;

error_phi_seivd_upper = quantile(chain(start_new:end,10:18),0.68) - quantile(chain(start_new:end,10:18),0.5);
error_phi_seivd_lower = quantile(chain(start_new:end,10:18),0.5) - quantile(chain(start_new:end,10:18),0.32);


error_tau_seivd = std(chain(start_new:end,19:27)) ;




%%


start_new_2 = 1;



load('./../../community/results/cba_seiv.mat');
chain_final(:,6:10) = exp(2.303.*chain_final(:,6:10));

beta_seiv = median(chain_final(3000:end,1:5));
phi_seiv = median(chain_final(3000:end,6:10));
tau_seiv = median(chain_final(3000:end,11:15));
r_seiv = median(chain_final(3000:end,16:18));


error_beta_seiv = std(chain_final(3000:end,1:5)) ;
error_phi_seiv = std(chain_final(3000:end,6:10)) ;
error_tau_seiv = std(chain_final(3000:end,11:15)) ;
error_r_seiv = std(chain_final(3000:end,16:18)) ;


load('./../../community/results/psa_seiv.mat');
chain_final(:,5:8) = exp(2.303.*chain_final(:,5:8));

beta_seiv(6:9) = median(chain_final(2000:end,1:4));
phi_seiv(6:9) = median(chain_final(2000:end,5:8));
tau_seiv(6:9) = median(chain_final(2000:end,9:12));
r_seiv(4:5) = median(chain_final(2000:end,13:14));

error_beta_seiv(6:9) = std(chain_final(2000:end,1:4));
error_phi_seiv(6:9) = std(chain_final(2000:end,5:8));
error_tau_seiv(6:9) = std(chain_final(2000:end,9:12));
error_r_seiv(4:5) = std(chain_final(2000:end,13:14));


%% one-steps

% one steps conventional 
beta_onestep_conventional = [0.94,91.2,25.7,0,10.4,44.3,319,54.2,357];
tau_onestep_conventional = [1,1.25,1.5,0,1,0.833,1.5,0.667,1.167];
%phi_onestep_conventional = [1.83e-7,1.88e-8,1.33e-7,0,9.9e-8,1.79e-7,7.6e-8,1.87e-7,1e-7];

% corrected for the indices.
phi_onestep_conventional = [1.88e-8,1.83e-7,1.33e-7,9.9e-8,1.79e-7,1.87e-7,7.6e-8,1e-7];

%one step bayesian (these are corrected for the indices)
% beta_onestep = [118.2,1.18,60.5,-200,8.246,99.2,93,437,413];
% tau_onestep = [1.89,1.749,2.19,0,1.9,1.47,1.42,2.14,1.97];
% phi_onestep = [ 1.45e-8,5.14e-8, 9.57e-8,0,1.227e-7,1.56e-7, 1.31e-7, 6.46e-8, 8.02e-8 ];
% 
% beta_onestep_error = [68.8,2.02,50,0,21.9,43.8,50,79.5,86];
% tau_onestep_error = [0.18,0.56,0.19,0,0.4,0.13,0.12,0.13,0.09];
% phi_onestep_error = [9.57e-9, 1.92e-8,1.39e-8,0,1.68e-8, 2.29e-8,1.99e-8,1.75e-8,1.7e-8];
% 
% r_onestep_conventional = [0.19,0.245,0.22,0.28,0.25];


% before replacing hs6_h100 with the proper file.
% beta_onestep = [   126.55        1.901       63.884     -10   36.55       75.092       87.229       351.51       324.15];
% tau_onestep = [   1.8781       1.9068       2.0352   -10    1.8012       1.4589       1.3709       2.1661       2.0193];
% phi_onestep = [ 1.776e-08   5.6449e-08   1.2386e-07  -10  3.506e-08   1.2251e-07   8.5023e-08   5.5815e-08   4.4429e-08];
% 
% 
% 
% beta_onestep_error = [    36.508      0.12939       48.102    0   45.319       30.193       42.213       17.777        71.87];
% tau_onestep_error = [   0.038254      0.26513      0.17063    0  0.27281      0.11118     0.075967     0.048039     0.078438];
% phi_onestep_error = [ 5.2604e-09   1.1261e-08   2.7217e-08  0 1.2592e-08   1.3211e-08    9.791e-09   8.1126e-09   7.4481e-09];
% 
% 
%  r_bayesian = [      0.2531      0.18702      0.23981      0.21843      0.28019      0.25393      0.27825      0.25793];
%  r_bayesian_std = [  0.013507     0.048525     0.045246     0.046714     0.042521     0.049102     0.040756     0.019837];




% new
beta_onestep = [   126.55        1.901       63.884    -10    36.55       75.092       87.229      431.35      324.15];
%beta_onestep = [   126.55        1.901       63.884    -10    36.55       75.092       87.229      357.25      324.15];
beta_onestep = [   126.55        1.901       63.884    -10    36.55       75.092       87.229       435.6      324.15];


tau_onestep = [   1.8781       1.9068       2.0352    -10   1.8012       1.4589       1.3709       2.22       2.0193];
phi_onestep = [ 1.776e-08   5.6449e-08   1.2386e-07  -10  3.506e-08   1.2251e-07   8.5023e-08   5.60e-08   4.4429e-08];



beta_onestep_error = [    36.508      0.12939       48.102    0   45.319       30.193       42.213        100        71.87];
%17.77
%beta_onestep_error = [    36.508      0.12939       48.102    0   45.319       30.193       42.213        17.77       71.87];

tau_onestep_error = [   0.038254      0.26513      0.17063    0  0.27281      0.11118     0.075967     0.10    0.078438];
phi_onestep_error = [ 5.2604e-09   1.1261e-08   2.7217e-08   0  1.2592e-08   1.3211e-08    9.791e-09    7.6718e-09   7.4481e-09];


 r_bayesian = [      0.2531      0.18702      0.23981      0.21843      0.28019      0.25393       0.27987      0.25793];
 r_bayesian_std = [  0.013507     0.048525     0.048     0.046714     0.042521     0.049102     0.040756     0.019837];




r_onestep = [r_bayesian(2), mean([r_bayesian(2),r_bayesian(3)]), r_bayesian(4), mean([r_bayesian(5),r_bayesian(7)]),mean([r_bayesian(6),r_bayesian(8)]) ];
r_onestep_error = [r_bayesian_std(2), sqrt(r_bayesian_std(2)^2 + r_bayesian_std(3)^2), r_bayesian_std(4),  sqrt(r_bayesian_std(5)^2 + r_bayesian_std(7)^2), sqrt(r_bayesian_std(6)^2 + r_bayesian_std(8)^2) ];







%%  one step vs community

gap = 0.1;
market_size_given = 10;
%color_green = [70/255,210/255,130/255];
color_green = [171,193,157]./255;

figure
subplot(2,2,1)
errorbar(x+gap,beta_onestep,beta_onestep_error,'MarkerSize',market_size_given,'MarkerFaceColor',[0.5,0.5,0.5],'Marker','o','Color',[0.5,0.5,0.5],'LineWidth',2,"LineStyle","none");
hold on;
errorbar(x-gap,beta,error_beta_seivd, 'o','MarkerSize',market_size_given,'MarkerFaceColor',color_green,'Color',color_green,'LineWidth',2);
%errorbar(x,beta_seiv,error_phi_seiv,'MarkerSize',market_size_given,'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','square','LineWidth',2, ...
%     'LineStyle','none',"Color",'k')

set(gca,'FontSize',20,'FontName','times');
xticks(1:9);
ylabel({'Burst sizes' ;'\beta (virions/cell)'});
ylim([0 600]);
yticks(0:100:1000);
xline(5.5,'--k',LineWidth=2.5);
%legend('Inferred parameters from scaled-up SEIVD model (Fig 4)','Inferred parameters from scaled-up SEIV model (Fig 3)','Inferred parameters from pairwise SEIV model (Fig 2)');
 %legend('Parameters inferred with scaled-up SEIVD community model','Parameters inferred with scaled-up SEIV community model','Parameters from pairwise SEIV model');
 %legend('Box','off');

subplot(2,2,2)
errorbar(x-gap,tau_onestep,tau_onestep_error,'MarkerSize',market_size_given,'MarkerFaceColor',[0.5,0.5,0.5],'Marker','o','Color',[0.5,0.5,0.5],'LineWidth',2,"LineStyle","none");
hold on;
errorbar(x+gap,tau,error_tau_seivd, 'o','MarkerSize',market_size_given,'MarkerFaceColor',color_green,'Color',color_green,'LineWidth',2);


%errorbar(x,tau_seiv,error_tau_seiv,'MarkerSize',market_size_given,'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','square','LineWidth',2, ...
%     'LineStyle','none',"Color",'k')
set(gca,'FontSize',20, 'FontName','times');
xticks(1:9);
ylabel({'Latent periods'; ' \tau (hr)'});
yticks([2 4 6 8 10]);
ylim([1 6]);
xline(5.5,'--k',LineWidth=2.5);
yticks(1:6)


% 
subplot(2,2,3)

errorbar(x-gap,phi_onestep, phi_onestep_error,'MarkerSize',market_size_given,'MarkerFaceColor',[0.5,0.5,0.5],'Marker','o','Color',[0.5,0.5,0.5],'LineWidth',2,"LineStyle","none");
hold on;
%errorbar(x+gap,phi,error_phi_seivd, 'o','MarkerSize',market_size_given,'MarkerFaceColor',color_green,'Color',color_green,'LineWidth',2);

errorbar(x+gap,phi,error_phi_seivd_lower,error_phi_seivd_upper, 'o','MarkerSize',market_size_given,'MarkerFaceColor',color_green,'Color',color_green,'LineWidth',2);


%errorbar(x,phi_seiv,error_phi_seiv,'MarkerSize',market_size_given,'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','square','LineWidth',2, ...
%     'LineStyle','none',"Color",'k')


set(gca,'FontSize',20,'FontName','times');
xticks(1:9);
set(gca, 'YScale', 'log');
ylabel({'Adsorption rates';'\phi (ml/hr)'});
%yticks([1e-10 1e-9 1e-8 1e-7]);
xline(5.5,'--k',LineWidth=2.5);
%ylim([1e-10 10e-7]);




subplot(2,2,4)
set(gca, 'YScale', 'linear')


errorbar((1:5)-gap,r_onestep,r_onestep_error,'MarkerSize',market_size_given,'MarkerFaceColor',[0.5,0.5,0.5],'Marker','o','Color',[0.5,0.5,0.5],'LineWidth',2,"LineStyle","none");
hold on;
errorbar((1:5)+gap,r,error_r_seivd, 'o','MarkerSize',market_size_given,'MarkerFaceColor',color_green,'Color',color_green,'LineWidth',2);




%errorbar(1:5,r_seiv,error_r_seiv,'MarkerSize',market_size_given,'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','square','LineWidth',2, ...
%   'LineStyle','none',"Color",'k')


set(gca,'FontSize',20,'FontName','times');
xticks(1:5);
xticklabels({'CBA 4','CBA 18','CBA 38','PSA H100','PSA 13-15'});
set(gca,'XTickLabelRotation',90)
ylabel({'Growth rates','r (cells/hr)'});
xline(3.5,'--k',LineWidth=2.5);
legend("Pairwise parameters","Community parameters")


