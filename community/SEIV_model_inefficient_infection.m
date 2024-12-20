

% just to plot the model with the inferred median parameters (approx)
% no dependancies -- these dependancies have driven me crazy!

clear;
clc;
%% Settings for running the model

addpath('./src/models/');

model = SEIV_diff_NE_ineffi_infection(5,5,200);
model.name = 'SEIV-noninfectious';
model.debris_inhib = 0;
model.debris_inhib2 = 0;
model.debris_inhib3 = 0;
model.debris_inhib4 = 0;
model.debris_inhib5 = 0;


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

pars_afterinf.epsilon = [1 1 1 1 1 1 1 1 1 1];

% not important -- the model class does not call these paramters, the
% plots should not change if these parameters are changed.
pars_afterinf.a = eye(5);
pars_afterinf.m = [  0.00081472
   0.00090579
   0.00012699
   0.00091338
   0.00063236];

pars_afterinf.q = [0.5
          0.5
          0.5
          0.5
          0.5];

pars_afterinf.prob = [    0
     0
     0
     0
     0];



%% simulate 

tvec = 0:0.1:15.75;
[t_after,S_after,V_after,D_after,I_after,E_after] =  simulate_ode(model,pars_afterinf,tvec,pars_afterinf.S0,pars_afterinf.V0); % mcmc parameter set


%%  load raw data and plot
load('./data/triplicate_data.mat');

clear ans;
color_ofthe_fit = [1 0 0]*0.5;
linewidth = 2;

hf4 = figure;
subplot(2,5,1)
errorbar(time/60,mean(1e3*host1'),std(1e3*host1'),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
ylabel({'Host density';'(cells/ml)'});
title('CBA 4','FontSize',18);
    plot(t_after,S_after(:,1),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,2)
errorbar(time/60,mean(1e3*host2'),std(1e3*host2'),'o','MarkerSize',8,  'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255] );hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 18','FontSize',18);
   
    plot(t_after,S_after(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,3)
errorbar(time/60,mean(1e3*host3'),std(1e3*host3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
title('CBA 38','FontSize',18);

    plot(t_after,S_after(:,3),'-','Color',color_ofthe_fit,'LineWidth',linewidth);



subplot(2,5,4)
errorbar(time/60,mean(1e3*host4'),std(1e3*host4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA H100','FontSize',18);
    
    plot(t_after,S_after(:,4),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,5)
errorbar(time/60,mean(1e3*host5'),std(1e3*host5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA 13-15','FontSize',18);
    
    plot(t_after,S_after(:,5),'-','Color',color_ofthe_fit,'LineWidth',linewidth);

%xlabel("Time (hours)");
%ylabel("Host density (cell/ml)");






subplot(2,5,6)
errorbar(time/60,mean(1e3*virus1'),std(1e3*virus1'),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e4 1e6 1e8 1e10]);
ylabel({'Phage density';'(virions/ml)'});
title('\phi18:2','FontSize',18);
    
    plot(t_after,V_after(:,1),'-','Color',color_ofthe_fit,'LineWidth',linewidth);

subplot(2,5,7)
errorbar(time/60,mean(1e3*virus2'),std(1e3*virus2'),'o','MarkerSize',8,  'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255] );hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e4 1e6 1e8 1e10]);
    title('\phi18:3','FontSize',18);
    
    plot(t_after,V_after(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,8)
errorbar(time/60,mean(1e3*virus3'),std(1e3*virus3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('\phi38:1','FontSize',18);
    
    plot(t_after,V_after(:,3),'-','Color',color_ofthe_fit,'LineWidth',linewidth);

subplot(2,5,9)
errorbar(time/60,mean(1e3*virus4'),std(1e3*virus4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('PSA HP1','FontSize',18);
    
    plot(t_after,V_after(:,4),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,10)
errorbar(time/60,mean(1e3*virus5'),std(1e3*virus5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],Color=[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('PSA HS6','FontSize',18);
    
    plot(t_after,V_after(:,5),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
    %legend('Data','95% confidence interval','Bayesian fit');
    %legend('Box','off');

han=axes(hf4,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
xlabel("Time (hours)");
