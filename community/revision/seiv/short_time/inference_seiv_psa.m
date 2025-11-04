%%% this one uses gaussian reparam trick -- 
%%% chooses a variance from 'sensitivity' of data and mean from point estimate of one experiment
%%% so definitely priors not reused.

clc;
clear all;

addpath(genpath('./../../..'));
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
mcmcpars = rmfield(mcmcpars,'epsilon');

%%

mcmcpars2.beta.subid = [1 2 3 4];
mcmcpars2.phi.subid = [1 2 3 4];
mcmcpars2.tau.subid = [1 2 3 4];
mcmcpars2.r.subid = [1 2];
mcmcpars2.r.log = 0;
mcmcpars2.beta.log = 0;
mcmcpars2.tau.log = 0;
mcmcpars2.phi.log = 1;



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

mcmcoptions.nsimu = 50000;
mcmcoptions.method  = 'dram';

transient_id = 30;
lambda = 0;


clear mcmcparam;
for i = 1:14
mcmcparam{1,i}{1,1} = 'not_named'; % feeling lazy :)
mcmcparam{1,i}{1,3} = -4;
mcmcparam{1,i}{1,4} = 4;

mcmcparam{1,i}{1,2} = 0;
mcmcparam{1,i}{1,5} = 0;
mcmcparam{1,i}{1,6} = 1;
end



theta_std = [100
    100
    50
    50
    0.2
    0.2
    0.2
    0.2
    0.2
    0.2
    0.2
    0.2
    0.1
    0.1];


mcmcmodel.sigma2 = 5; % (initial) error variance from residuals of the lsq fit
mcmcmodel.S20 = mcmcmodel.sigma2;
mcmcmodel.N0 = 10;
mcmcoptions.updatesigma = 1;

limit = 0.25;
mcmcmodel.ssfun = @(theta,data) error_psa_reparam(theta,data,pars2,model,mcmcpars2,limit,theta_initial,theta_std); 
[mcmcresults, chain, s2chain]= mcmcrun(mcmcmodel,data,mcmcparam,mcmcoptions);




chain_final = repmat(theta_initial',mcmcoptions.nsimu,1) + chain.*repmat(theta_std',mcmcoptions.nsimu,1);


%% save
save("psa_seiv_2.mat");

%%

% figure(2)
% 
% subplot(2,2,1)
% errorbar(time/60,mean(1e3*host4'),std(1e3*host4'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],Color = [70/255,130/255,180/255]);hold on;
% set(gca, 'YScale', 'log');set(gca,'FontSize',20)
% set(gca,'fontname','times')  % Set it to times
% ylim([1e5 1e8]);
%     xlim([0 16]);
%   xticks([0 2 4 6 8 10 12 14 16]);
%   axis('square');
%     yticks([1e5 1e6 1e7 1e8]);
%     title('PSA H100','FontSize',18);
%     
%   
% 
% subplot(2,2,2)
% errorbar(time/60,mean(1e3*host5'),std(1e3*host5'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],Color = [70/255,130/255,180/255]);hold on;
% set(gca, 'YScale', 'log');set(gca,'FontSize',20);
% set(gca,'fontname','times')  % Set it to times
% ylim([1e5 1e8]);
%     xlim([0 16]);
%     xticks([0 2 4 6 8 10 12 14 16]);
%     axis('square');
%     yticks([1e5 1e6 1e7 1e8]);
%     title('PSA 13-15','FontSize',18);
%     
%   
% 
% subplot(2,2,3)
% errorbar(time/60,mean(1e3*virus4'),std(1e3*virus4'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],Color = [70/255,130/255,180/255]);hold on;
% set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
% ylim([1e4 1e11]);
%     xlim([0 16]);
%    xticks([0 2 4 6 8 10 12 14 16]);
%     set(gca,'FontSize',20);
%     axis('square');
%    yticks([1e4 1e6 1e8 1e10]);
%    title('PSA HP1','FontSize',18);
%    
% 
% 
% subplot(2,2,4)
% errorbar(time/60,mean(1e3*virus5'),std(1e3*virus5'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],Color = [70/255,130/255,180/255]);hold on;
% set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
% ylim([1e4 1e11]);
%     xlim([0 16]);
%  xticks([0 2 4 6 8 10 12 14 16]);
%     set(gca,'FontSize',20);
%     axis('square');
%   yticks([1e4 1e6 1e8 1e10]);
%   title('PSA HS6','FontSize',18);
%     
% 
% for i = 3000:10:5000
% 
%     i
%   if error_psa_short(chain_final(i,:),pars2,data,model,mcmcpars2,0.5) ~= Inf
% pars2 = update_pars(pars2,chain_final(i,:),mcmcpars2);
% pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);
% [t2,S_median,V_median,D_median,I_median,E_median] =  simulate_ode(model,pars2,tvec,pars2.S0,pars2.V0); % mcmc parameter set
% 
% subplot(2,2,1)
% patchline(t2,S_median(:,1),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);
% 
% 
% subplot(2,2,2)
% patchline(t2,S_median(:,2),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);
% 
% 
% subplot(2,2,3)
%   patchline(t2,V_median(:,1),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);
% 
% 
% subplot(2,2,4)
% patchline(t2,V_median(:,2),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);
%  
%   end
% 
% end
% 
% 
