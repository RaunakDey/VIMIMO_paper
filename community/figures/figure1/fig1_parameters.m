clear all;
clc;

x=1:8;




%% one-steps

% one steps conventional 
% beta_onestep_conventional = [0.94,91.2,25.7,10.4,44.3,319,54.2,357];
% tau_onestep_conventional = [1,1.25,1.5,1,0.833,1.5,0.667,1.167];
% phi_onestep_conventional = [1.83e-7,1.88e-8,1.33e-7,9.9e-8,1.79e-7,7.6e-8,1.87e-7,1e-7];


beta_onestep_conventional = [91.2,0.94,27.29,10.5,54.9,54.2,205.8,276.7];
tau_onestep_conventional = [1.25,1,1.41,1,0.75,0.667,1.33,1.25];
phi_onestep_conventional = [1.88e-8,1.83e-7,1.33e-7,1.65e-8,9.9e-8,1.79e-7,1.87e-7,7.6e-8,1e-7];

beta_onestep_conventional_error = [32.91,0.15,6.28,3.18,21.8,15.74,92.7,110.4];
phi_one_step_conventional_error = [3.4e-09     3.99e-08     2.66e-08     4.96e-10     1.37e-08     3.33e-08     3.79e-08     8.64e-09     3.55e-09]';


   r_onestep_conventional_error =   [0.019093     0.021355     0.030257     0.052883     0.044104]';

%one step bayesian
% beta_onestep = [1.81,118.2,60.5,8.246,99.2,437,93,413];
% tau_onestep = [1.749,1.89,2.19,1.9,1.47,2.14,1.42,1.97];
% phi_onestep = [5.14e-8,1.45e-8,9.57e-8,1.227e-7,1.56e-7, 6.46e-8,1.31e-7,8.02e-8 ];
% 
% beta_onestep_error = [2.02,68.8,50,21.9,43.8,79.5,50,86];
% tau_onestep_error = [0.56,0.18,0.19,0.4,0.13,0.13,0.12,0.09] ;
% phi_onestep_error = [1.92e-8,9.57e-9,1.39e-8,1.68e-8, 2.29e-8,1.75e-8,1.99e-8,1.7e-8];

beta_onestep = [118.2,1.18,60.5,8.246,99.2,93,437,413];
tau_onestep = [1.89,1.749,2.19,1.9,1.47,1.42,2.14,1.97];
phi_onestep = [ 1.45e-8,5.14e-8, 9.57e-8,1.227e-7,1.56e-7, 1.31e-7, 6.46e-8, 8.02e-8 ];

beta_onestep_error = [68.8,2.02,50,21.9,43.8,50,79.5,86];
tau_onestep_error = [0.18,0.56,0.19,0.4,0.13,0.12,0.13,0.09];
phi_onestep_error = [9.57e-9, 1.92e-8,1.39e-8,1.68e-8, 2.29e-8,1.99e-8,1.75e-8,1.7e-8];







r_onestep_conventional = [0.19,0.245,0.22,0.28,0.25];

boxes_onestep = [77,74,85,55,83,98,87,109];
boxes_onestep_error = [41,39,37,40,34,35,35,31];

%% figures -- one step

%color_def = [70/255,130/255,180/255];
%color_def = [70/255,210/255,130/255];

color_def2 = [0.9290 0.6940 0.1250];
color_def = [111,193,157]./255;


figure(2)
subplot(1,5,1)
errorbar(x-0.1,beta_onestep_conventional,beta_onestep_conventional_error,'^','MarkerSize',14,'MarkerFaceColor',color_def2,'Color',color_def2);
hold on;
errorbar(x+0.1,beta_onestep,beta_onestep_error,'bo','MarkerSize',14,'MarkerFaceColor',color_def,'Color',color_def);
title({'Burst size', '\beta (virions/cell)'});
ylim([1 600]);
set(gca,'FontSize',18);
xticks(1:8)
xticklabels([1,2,3,5,6,7,8,9]);
yticks(0:100:600);
ylim([-1 600])
xlim([0 9]);


subplot(1,5,2)
plot(x-0.1,tau_onestep_conventional,'^','MarkerSize',14,'MarkerFaceColor',color_def2,'Color',color_def2);
hold on;
errorbar(x+0.1,tau_onestep,tau_onestep_error,'bo','MarkerSize',14,'MarkerFaceColor',color_def,'Color',color_def);
title({'Latent period',' \tau (hr)'});
ylim([0.2 3]);
set(gca,'FontSize',18);
xticks(1:8);
xticklabels([1,2,3,5,6,7,8,9]);
xlim([0 9]);


subplot(1,5,3)
errorbar([1 2 3 4 5 6 7 8 9]-0.1,phi_onestep_conventional,phi_one_step_conventional_error,'^','MarkerSize',14,'MarkerFaceColor',color_def2,'Color',color_def2);
hold on;
errorbar([1 2 3 5 6 7 8 9]+0.1,phi_onestep,phi_onestep_error,'bo','MarkerSize',14,'MarkerFaceColor',color_def,'Color',color_def);
title({'Adsorption rate','\phi (ml/hr)'})
ylim([1e-10 2.5e-7]);
set(gca,'FontSize',18);
xticks(1:9)
xticklabels([1,2,3,4,5,6,7,8,9]);
xlim([0 10]);
%legend('Conventional','Bayesian');
%legend('boxoff');



cv_mean_onestep  = 1./sqrt(boxes_onestep);
cv_onestep_error = 0.5* (boxes_onestep).^(-1.5).*boxes_onestep_error;
%cv_onestep_error  = 1./sqrt(boxes_onestep + boxes_onestep_error) - 1./sqrt(boxes_onestep + boxes_onestep_error)

subplot(1,5,4)
errorbar(x,cv_mean_onestep,cv_onestep_error,'ko','MarkerSize',14,'MarkerFaceColor',color_def,'Color',color_def);
title({'Coefficient of',' variation of \tau (CV)'});
set(gca,'FontSize',18);
xticks(1:8);
xticklabels([1,2,3,5,6,7,8,9]);
ylim([0 0.25])
%yticks(0:0.1:0.3);
xlim([0 9]);


subplot(1,5,5)
errorbar(1:5,r_onestep_conventional,r_onestep_conventional_error,'^','MarkerSize',14,'MarkerFaceColor',color_def2,'Color',color_def2);
title({'Growth rate','r (cells/hr)'});
ylim([0.15 0.35]);
set(gca,'FontSize',18);
xlim([0.5 5.5])
xticks(1:5);
xticklabels({'CBA 4','CBA 18','CBA 38','PSA H100','PSA 13-15'});
set(gca,'XTickLabelRotation',90);




