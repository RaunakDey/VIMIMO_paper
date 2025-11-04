clc;
clear all;

addpath(genpath('./../..'))
load('./../../combined_posteriors.mat');

load('./../../../data/triplicate_data.mat');
load('./../../seiv-data.mat');

linewidth = 2;

color_ofthe_fit = [1 0 0]*0.5;
color_ofthe_fit = [0 0 0];
color_ofthe_fill = [0.95 0 0];
transparency = 0.25;
tvec = 0:0.01:15.75; % for better viz






clear pars2
%% old parameters


pars2.epsilon = ones(1,4); %hard coding off multiplicative correction bias.

pars2.phi =    1.0e-07 *[   0.63   0.31;
           0.63   0.31];

pars2.beta = [700  50;
        700  50];

pars2.tau = [  10    10;
          10  10];

pars2.eta = zeros(2,2);
pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);

pars2.r = [ 0.7   0.5]';

pars2.NE = [83    70;
      87   70];

pars2.S0 = 1e6*[ 6.2050
    7.7533];

pars2.V0 = 1e7*[ 0.0110
    1.1510];

pars2.NH = 2;
pars2.NV = 2;
pars2.M = ones(2,2);

max_NE = round(max(max(pars2.NE)));
model = SEIVD_diff_NE_diff_debris_abs(2,2,max_NE);
model.host_growth = 0;
model.viral_decay = 0;
model.viral_adsorb = 0;
model.lysis_reset = 0;
model.debris_inhib = 0;
model.debris_inhib2 = 0;


%% find min -- point estimate

mcmcpars = rmfield(mcmcpars,'epsilon');

theta_initial = [pars2.beta(1,1)
    pars2.beta(2,1)
    pars2.beta(1,2)
    pars2.beta(2,2)
    log(pars2.phi(1,1))/log(10)
    log(pars2.phi(2,1))/log(10)
    log(pars2.phi(1,2))/log(10)
    log(pars2.phi(2,2))/log(10)
    pars2.tau(1,1)
    pars2.tau(2,1)
    pars2.tau(1,2)
    pars2.tau(2,2)
    pars2.r(1)
    pars2.r(2)
    ];

mcmcpars2.beta.subid = [1 2 3 4];
mcmcpars2.phi.subid = [1 2 3 4];
mcmcpars2.tau.subid = [1 2 3 4];
mcmcpars2.r.subid = [1 2];
mcmcpars2.r.log = 0;
mcmcpars2.beta.log = 0;
mcmcpars2.tau.log = 0;
mcmcpars2.phi.log = 1;
%%

theta_initial= [    960.0000
  960.0000
  190.0000
   10.0000
   -8.2007
   -8.0007
   -8.5086
   -9.1086
    8.4000
    9.8000
    8.2000
   10.0000
    0.1500
    0.0500];



% theta_optimized = [    1110
%          1080
%           280
%           -60
%       -7.0007
%       -8.6007
%       -8.1086
%       -9.3086
%           5.4
%           6.8
%           5.2
%             7
%          0.45
%           0.3];


% theta_optimized =[ 1110
%          1070
%           300
%             0
%       -7.0007
%       -8.6007
%       -8.1086
%       -9.5086
%           5.4
%           6.8
%           5.2
%             7
%           0.5
%          0.25];


% theta_optimized = [ 700
%   700
%   50
%   50
%   -7.2
%    -7.2
%    -7.7
%    -7.7
%     10
%     10
%     3
%     3
%     0.7
%     0.5];



theta_initial = [ 660.0000
  600.0000
  70.0
  70.0
   -7.2
   -7.2
   -7.8
   -8.2
    4.2875
    4.3115
    3.7246
   3.9793
    0.7
    0.5];

     theta_initial =   [ 550
          750
           40
           70
         -7.2
         -7.2
         -7.8
         -8.2
       3.0875
       5.9115
       3.5246
       5.1793
         0.75
         0.68];


     theta_initial =[400
          900
           40
           70
         -7.2
         -7.2
         -7.8
         -8.2
       3.2875
       5.1115
       3.3246
       5.7793
         0.75
         0.68];

theta_sd = [10
    10
    10
    10
    0.2
    0.2
    0.2
    0.2
    0.2
    0.2
    0.2
    0.2
    0.05
    0.05];





pars2 = update_pars(pars2,theta_initial,mcmcpars2);

pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);

error_psa_short(theta_initial,pars2,data,model,mcmcpars2,0.25)

theta_optimized = search_minimum_psa(theta_initial,theta_sd,data,model,pars2,mcmcpars2,0.25,15,1);

error_psa_short(theta_optimized,pars2,data,model,mcmcpars2,0.25)

pars2 = update_pars(pars2,theta_optimized,mcmcpars2);

pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);

[t2,S_median,V_median,D_median,I_median,E_median] =  simulate_ode(model,pars2,tvec,pars2.S0,pars2.V0); % mcmc parameter set

%%
hf4 = figure;


subplot(2,2,1)
errorbar(time/60,mean(1e3*host4'),std(1e3*host4'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA H100','FontSize',18);
    
    plot(t2,S_median(:,1),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,2,2)
errorbar(time/60,mean(1e3*host5'),std(1e3*host5'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA 13-15','FontSize',18);
    
    plot(t2,S_median(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth);

%xlabel("Time (hours)");
%ylabel("Host density (cell/ml)");







subplot(2,2,3)
errorbar(time/60,mean(1e3*virus4'),std(1e3*virus4'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('PSA HP1','FontSize',18);
    
    plot(t2,V_median(:,1),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,2,4)
errorbar(time/60,mean(1e3*virus5'),std(1e3*virus5'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('PSA HS6','FontSize',18);
    
    plot(t2,V_median(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
    %legend('Data','95% confidence interval','Bayesian fit');
    %legend('Box','off');

han=axes(hf4,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
xlabel("Time (hours)");

