clc;
clear all;

addpath(genpath('./../..'))
load('./../combined_posteriors.mat');

load('./../../data/triplicate_data.mat');
load('./../seiv-data.mat');

linewidth = 2;

color_ofthe_fit = [1 0 0]*0.5;
color_ofthe_fit = [0 0 0];
color_ofthe_fill = [0.95 0 0];
transparency = 0.25;
tvec = 0:0.01:15.75; % for better viz







%% old parameters


pars2.epsilon = ones(1,10); %hard coding off multiplicative correction bias.

pars2.phi =    1.0e-07 *[         0    0.0015         0         0        0;
    0.0876    0.0906    0.0313         0         0;
         0         0    0.0063         0         0;
         0         0         0    0.0218    0.1835;
         0         0         0    0.1193    0.0248];

pars2.beta = [0  359.7693         0         0         0;
  126.4146  436.8215  313.7494         0         0;
         0         0  107.3469         0         0;
         0         0         0  192.1766  112.9907;
         0         0         0  426.0895   21.0476];

pars2.tau = [ 0    6.7845         0         0         0;
    4.0770    4.6333    4.4117         0         0;
         0         0    3.5320         0         0;
         0         0         0    8.1820    5.4737;
         0         0         0    1.8312   18.9541];

pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);

pars2.r = [0.1   0.1726    0.0107    0.4997    0.5874]';

pars2.NE = [0    77     0     0     0;
    74    85    71     0     0;
     0     0    55     0     0;
     0     0     0    83    98;
     0     0     0    87   109];


%% new parameters

pars2.epsilon = ones(1,10); %hard coding off multiplicative correction bias.

pars2.phi =    1.0e-07 *[  0    0.65         0         0        0;
                           0.176    0.806    0.1         0         0;
                             0         0    0.163         0         0;
                             0         0         0    0.518    0.1835;
                             0         0         0    0.5193    0.0248];

pars2.beta = [0    4       0         0         0;
             500  200    310         0         0;
              0      0     107         0         0;
              0         0         0  192 112;
              0         0         0  420   21];

pars2.tau = [ 0   5        0         0         0;
    2.7    5    6.4117         0         0;
         0         0    3.5320         0         0;
         0         0         0    8.1820    5.4737;
         0         0         0    1.8312   18.9541];

pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);

pars2.r = [0.1894    0.1726    0.307    0.65    0.5874]';

pars2.NE = [0    30     0     0     0;
    10    80    71     0     0;
     0     0    55     0     0;
     0     0     0    50   50;
     0     0     0    50   50];



max_NE = round(max(max(pars2.NE)));
model = SEIVD_diff_NE_diff_debris_abs(5,5,max_NE);
model.host_growth = 0;
model.viral_decay = 0;
model.viral_adsorb = 0;
model.lysis_reset = 0;
model.debris_inhib = 0;
model.debris_inhib2 = 0;
model.debris_inhib3 = 0;
model.debris_inhib4 = 0;
model.debris_inhib5 = 0;

%% find min -- point estimate

mcmcpars = rmfield(mcmcpars,'epsilon');

theta_initial(1:5) = pars2.r;
theta_initial(6:14) = pars2.beta(pars2.beta>0);
theta_initial(15:23) = log(pars2.phi(pars2.phi>0))/log(10);
theta_initial(24:32) = pars2.tau(pars2.tau>0);




theta_optimized =[0.2573
    0.1926
    0.2499
    0.65
    0.5864
  500.0000
    100.0001
  400
  310.0000
  107.0000
  191.9998
  419.9993
  111.9973
   20.9994
   -7.6604
   -6.8128
   -7.3191
   -7.9
   -7.7765
   -7.0942
   -8.0788
   -8.0752
   -8.8004
    2.0421
    8.4164
    6.3850
    6.4149
    3.5301
    4.2875
    4.3115
    3.7246
   3.9793]';



pars2 = update_pars(pars2,theta_optimized,mcmcpars);
pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);

[t2,S_median,V_median,D_median,I_median,E_median] =  simulate_ode(model,pars2,tvec,pars2.S0,pars2.V0); % mcmc parameter set
figure(1)
subplot(2,2,1)
plot(t2,sum(I_median,2)) %-- time series of sum of all infected cells.
ylabel('Infected cells (all hosts)');

subplot(2,2,2)
plot(t2,sum(I_median,2)+sum(E_median,2) + sum(S_median,2) )
ylabel('All hosts')

subplot(2,2,3)
plot(t2,D_median);
ylabel('Num dead cells (All hosts)')

%error_from_pars(pars2,data,model)

%%
hf4 = figure;
subplot(2,5,1)
errorbar(time/60,mean(1e3*host1'),std(1e3*host1'),'o','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
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
    plot(t2,S_median(:,1),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,2)
errorbar(time/60,mean(1e3*host2'),std(1e3*host2'),'o','MarkerSize',4,  'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255] );hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 18','FontSize',18);
   
    plot(t2,S_median(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,3)
errorbar(time/60,mean(1e3*host3'),std(1e3*host3'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
title('CBA 38','FontSize',18);

    plot(t2,S_median(:,3),'-','Color',color_ofthe_fit,'LineWidth',linewidth);



subplot(2,5,4)
errorbar(time/60,mean(1e3*host4'),std(1e3*host4'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA H100','FontSize',18);
    
    plot(t2,S_median(:,4),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,5)
errorbar(time/60,mean(1e3*host5'),std(1e3*host5'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA 13-15','FontSize',18);
    
    plot(t2,S_median(:,5),'-','Color',color_ofthe_fit,'LineWidth',linewidth);

%xlabel("Time (hours)");
%ylabel("Host density (cell/ml)");






subplot(2,5,6)
errorbar(time/60,mean(1e3*virus1'),std(1e3*virus1'),'o','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e4 1e6 1e8 1e10]);
ylabel({'Phage density';'(virions/ml)'});
title('\phi18:2','FontSize',18);
    
    plot(t2,V_median(:,1),'-','Color',color_ofthe_fit,'LineWidth',linewidth);

subplot(2,5,7)
errorbar(time/60,mean(1e3*virus2'),std(1e3*virus2'),'o','MarkerSize',4,  'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255] );hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e4 1e6 1e8 1e10]);
    title('\phi18:3','FontSize',18);
    
    plot(t2,V_median(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,8)
errorbar(time/60,mean(1e3*virus3'),std(1e3*virus3'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('\phi38:1','FontSize',18);
    
    plot(t2,V_median(:,3),'-','Color',color_ofthe_fit,'LineWidth',linewidth);

subplot(2,5,9)
errorbar(time/60,mean(1e3*virus4'),std(1e3*virus4'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('PSA HP1','FontSize',18);
    
    plot(t2,V_median(:,4),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,10)
errorbar(time/60,mean(1e3*virus5'),std(1e3*virus5'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('PSA HS6','FontSize',18);
    
    plot(t2,V_median(:,5),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
    %legend('Data','95% confidence interval','Bayesian fit');
    %legend('Box','off');

han=axes(hf4,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
xlabel("Time (hours)");

