%%% this one uses gaussian reparam trick -- 
%%% chooses a variance from 'sensitivity' of data and mean from point estimate of one experiment
%%% so definitely priors not reused.

clc;
clear all;

addpath(genpath('./../../..'))
addpath('./../../../src/models/');
addpath('./../../../mcmcstat/');
load('./../../combined_posteriors.mat');

load('./../../../data/triplicate_data.mat');

linewidth = 2;

color_ofthe_fit = [0.6 0.6 0.6];
color_ofthe_fill = [0.8 0 0];
transparency = 0.1;
tvec = 0:0.01:15.75; % for better viz


clear pars2
%% old parameters


pars2.epsilon = ones(1,6); %hard coding off multiplicative correction bias.



 pars2.phi = [0   10^(-7.2487)  0;
10^(-7.9454)  10^(-7.5) 10^(-8.3433);
0  0  10^(-7.6)];



pars2.beta = [0 100 0 ;
              500  400 310;
                0 0 106];


pars2.tau = [0 10 0;
            2 6 6;
            0 0  4];

pars2.eta = zeros(3,3);
pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);

pars2.r = [0.123 0.1873 0.24]';

pars2.NE = [0 30 0;
    10 80 71;
    0 0 60];

pars2.S0 =   1.0e+06 *[2.5111
    5.6423
    3.0257];

pars2.V0 =    1.0e+07 * [ 0.0429
    0.0287
    0.0528];

pars2.NH = 3;
pars2.NV = 3;
pars2.M = [0 1 0;
    1 1 1;
    0 0 1];



load('./cba_seiv.mat');
max_NE = round(max(max(pars2.NE)));
model = SEIVD_diff_NE_diff_debris_abs(3,3,max_NE);
model.host_growth = 0;
model.viral_decay = 0;
model.viral_adsorb = 0;
model.lysis_reset = 0;
model.debris_inhib = 0;
model.debris_inhib2 = 0;
model.debris_inhib3 = 0;

%mcmcpars = rmfield(mcmcpars,'epsilon');

%%

mcmcpars2.beta.subid = [2 4 5 8 9];
mcmcpars2.phi.subid = [2 4 5 8 9];
mcmcpars2.tau.subid = [2 4 5 8 9];
mcmcpars2.r.subid = [1 2 3];
mcmcpars2.r.log = 0;
mcmcpars2.beta.log = 0;
mcmcpars2.tau.log = 0;
mcmcpars2.phi.log = 1;


theta_initial = [500
   100
  400
  310
  106
   -7.9454
   -7.2487
   -7.5
   -8.3433
   -7.6
    2
    10
    6
  6
   4
    0.123
    0.1873
    0.24];


mcmcoptions.nsimu = 50000;
mcmcoptions.method  = 'dram';

transient_id = 30;
lambda = 0;


clear mcmcparam;
for i = 1:length(theta_initial)
mcmcparam{1,i}{1,1} = 'not_named'; % feeling lazy :)
mcmcparam{1,i}{1,3} = -4;
mcmcparam{1,i}{1,4} = 4;

mcmcparam{1,i}{1,2} = 0;
mcmcparam{1,i}{1,5} = 0;
mcmcparam{1,i}{1,6} = 1;
end



theta_std = [10
    10
    10
    10
    10
    0.2
    0.2
    0.2
    0.2
    0.2
    0.1
    0.1
    0.1
    0.1
    0.1
    0.05
    0.05
    0.05];



mcmcmodel.sigma2 = 25; % (initial) error variance from residuals of the lsq fit
mcmcmodel.S20 = mcmcmodel.sigma2;
mcmcmodel.N0 = 10;
mcmcoptions.updatesigma = 1;

limit = 0.4;
mcmcmodel.ssfun = @(theta,data) error_cba_reparam(theta,data,pars2,model,mcmcpars2,limit,theta_initial,theta_std); 
[mcmcresults, chain, s2chain]= mcmcrun(mcmcmodel,data,mcmcparam,mcmcoptions);




chain_final = repmat(theta_initial',mcmcoptions.nsimu,1) + chain.*repmat(theta_std',mcmcoptions.nsimu,1);


%% save
save("cba_seiv_2.mat");


%%

figure(2)


subplot(2,3,1)
errorbar(time/60,mean(1e3*host1'),std(1e3*host1'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255], Color=[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 4','FontSize',18);
    
  

subplot(2,3,2)
errorbar(time/60,mean(1e3*host2'),std(1e3*host2'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255], Color=[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 18','FontSize',18);
    


subplot(2,3,3)
errorbar(time/60,mean(1e3*host3'),std(1e3*host3'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255], Color=[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 38','FontSize',18);
    

  

subplot(2,3,4)
errorbar(time/60,mean(1e3*virus1'),std(1e3*virus1'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255], Color=[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('\phi18:1','FontSize',18);
   


subplot(2,3,5)
errorbar(time/60,mean(1e3*virus2'),std(1e3*virus2'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255], Color=[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('\phi18:3','FontSize',18);
    

subplot(2,3,6)
errorbar(time/60,mean(1e3*virus3'),std(1e3*virus3'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255], Color=[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('\phi38:1','FontSize',18);





for i = 3000:100:10000

    i
pars2 = update_pars(pars2,chain_final(i,:),mcmcpars2);
pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);
[t2,S_median,V_median,D_median,I_median,E_median] =  simulate_ode(model,pars2,tvec,pars2.S0,pars2.V0); % mcmc parameter set

subplot(2,3,1)
patchline(t2,S_median(:,1),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);


subplot(2,3,2)
patchline(t2,S_median(:,2),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);


subplot(2,3,3)
  patchline(t2,S_median(:,3),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);


subplot(2,3,4)
patchline(t2,V_median(:,1),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);



subplot(2,3,5)
  patchline(t2,V_median(:,2),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);


subplot(2,3,6)
patchline(t2,V_median(:,3),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);
 


end


