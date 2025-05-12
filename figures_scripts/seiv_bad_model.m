
load('./../community/results/cba_seiv_2.mat');
%%
color_ofthe_fit = [0.8 0.8 0.8];
color_ofthe_fit = [0 0 0];
tvec = 0:0.1:15.75;


red1 = [255,186,186]/255;
red2 = [255,123,123]/255;
red3 = [255,82,82]/255;
red4 = [255,0,0]/255;
red5 = [167,0,0]/255;


blue1 = [179,205,224]/255;
blue2 = [100,151,177]/255;
blue3 = [0,91,150]/255;
blue4 = [3,10,208]/255;
blue5 = [0,0,75]/255;


%%
figure(2)
subplot(2,5,1)
errorbar(time/60,mean(1e3*host1'),std(1e3*host1'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue1, Color=blue1);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 4','FontSize',18);
    
  

subplot(2,5,2)
errorbar(time/60,mean(1e3*host2'),std(1e3*host2'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue2, Color=blue2);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 18','FontSize',18);
    


subplot(2,5,3)
errorbar(time/60,mean(1e3*host3'),std(1e3*host3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue3, Color=blue3);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 38','FontSize',18);
    

  

subplot(2,5,6)
errorbar(time/60,mean(1e3*virus1'),std(1e3*virus1'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red1, Color=red1);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('\phi18:1','FontSize',18);
   


subplot(2,5,7)
errorbar(time/60,mean(1e3*virus2'),std(1e3*virus2'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red2, Color=red2);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('\phi18:3','FontSize',18);
    

subplot(2,5,8)
errorbar(time/60,mean(1e3*virus3'),std(1e3*virus3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red3, Color=red3);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('\phi38:1','FontSize',18);



%% draw samples

%upper_limit = [1000 1000 1000 1000 1000 -7 -7 -7 -7 -7 10 10 10 10 10 0.3 0.3 0.3];
%lower_limit = [10 10 10 10 10 -9 -9 -9 -9 -9 0.5 0.5 0.5 0.5 0.5 0.1 0.1 0.1];

upper_limit = [500 500 500 500 500 -7 -7 -7 -7 -7 5 5 5 5 5 0.3 0.3 0.3];
lower_limit = [10 10 10 10 10 -9 -9 -9 -9 -9 0.5 0.5 0.5 0.5 0.5 0.1 0.1 0.1];

draws = 500;

sampled_params = lhs_sampling(lower_limit,upper_limit, draws);


%%
%transparency = 0.9;
transparency = 0.1;
linewidth = 6;
skip = 1;

error_cba = zeros(1,draws);

for i = 1:draws

fprintf("Simulation %i \n",i)

pars2 = update_pars(pars2,sampled_params(i,:),mcmcpars2);
pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);
[t2,S_median,V_median,D_median,I_median,E_median] =  simulate_ode(model,pars2,tvec,pars2.S0,pars2.V0); % mcmc parameter set

error_cba(i) = error_cba_short(sampled_params(i,:),pars2,data,model,mcmcpars2,1);

figure(2)
subplot(2,5,1)
patchline(t2,S_median(:,1),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);


subplot(2,5,2)
patchline(t2,S_median(:,2),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);


subplot(2,5,3)
  patchline(t2,S_median(:,3),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);


subplot(2,5,6)
patchline(t2,V_median(:,1),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);



subplot(2,5,7)
  patchline(t2,V_median(:,2),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);


subplot(2,5,8)
patchline(t2,V_median(:,3),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);
 


end

%%

figure(2)
subplot(2,5,1)
errorbar(time/60,mean(1e3*host1'),std(1e3*host1'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue1, Color=blue1);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 4','FontSize',18);
    
  

subplot(2,5,2)
errorbar(time/60,mean(1e3*host2'),std(1e3*host2'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue2, Color=blue2);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 18','FontSize',18);
    


subplot(2,5,3)
errorbar(time/60,mean(1e3*host3'),std(1e3*host3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue3, Color=blue3);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 38','FontSize',18);
    

  

subplot(2,5,6)
errorbar(time/60,mean(1e3*virus1'),std(1e3*virus1'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red1, Color=red1);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('\phi18:1','FontSize',18);
   


subplot(2,5,7)
errorbar(time/60,mean(1e3*virus2'),std(1e3*virus2'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red2, Color=red2);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('\phi18:3','FontSize',18);
    

subplot(2,5,8)
errorbar(time/60,mean(1e3*virus3'),std(1e3*virus3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red3, Color=red3);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('\phi38:1','FontSize',18);


  %% good parameters

  pars3 = pars2;
  pars3.phi = 1.0e-07 *[0    0.2070         0
    0.1457    0.0107    0.0824
         0         0    0.0256];

  pars3.beta = [0  353.3034         0
  246.5709  165.5261  223.0705
         0         0  169.5272];

  pars3.r = [0.2209
    0.1322
    0.1732];

pars3.tau = [0    2.7869         0
    3.5886    1.7665    4.6805
         0         0    1.960];

pars3.eta= [0    0.3588         0
    0.2787    0.5661    0.2137
         0         0    0.5101];

pars3.NE = [0    30     0
    10    80    71
     0     0    60];

%%






%%

load("./../community/results/psa_seiv_2.mat");
tvec = 0:0.1:15.75;
color_ofthe_fit = [0.8 0.8 0.8];
linewidth = 6;

figure(2)
subplot(2,5,4)
errorbar(time/60,mean(1e3*host4'),std(1e3*host4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue4,Color = blue4);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA H100','FontSize',18);



subplot(2,5,5)
errorbar(time/60,mean(1e3*host5'),std(1e3*host5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue5,Color = blue5);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA 13-15','FontSize',18);



subplot(2,5,9)
errorbar(time/60,mean(1e3*virus4'),std(1e3*virus4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red4,Color = red4);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('PSA HP1','FontSize',18);



subplot(2,5,10)
errorbar(time/60,mean(1e3*virus5'),std(1e3*virus5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red5,Color = red5);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('PSA HS6','FontSize',18);


  %%

  upper_limit = [700 700 700 700 -7 -7 -7 -7 10 10 10 10 0.7 0.7];
lower_limit = [10 10 10 10 -9 -9 -9 -9 0.5 0.5 0.5 0.5 0.3 0.3];
draws = 500;

sampled_params = lhs_sampling(lower_limit,upper_limit, draws);

error_psa = [];
%%
transparency = 0.1;

tvec = 0:0.1:15.75;
for i = 1:draws

    i
    error_psa(i) = error_psa_short(sampled_params(i,:),pars2,data,model,mcmcpars2,0.43);

  if error_psa_short(sampled_params(i,:),pars2,data,model,mcmcpars2,0.43) ~= Inf
pars2 = update_pars(pars2,sampled_params(i,:),mcmcpars2);
pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);
[t2,S_median,V_median,D_median,I_median,E_median] =  simulate_ode(model,pars2,tvec,pars2.S0,pars2.V0); % mcmc parameter set

subplot(2,5,4)
patchline(t2,S_median(:,1),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);
hold on;

subplot(2,5,5)
patchline(t2,S_median(:,2),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);
hold on;


subplot(2,5,9)
  patchline(t2,V_median(:,1),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);
hold on;


subplot(2,5,10)
patchline(t2,V_median(:,2),'edgecolor',color_ofthe_fit,'edgealpha',transparency,'LineWidth',linewidth);
 hold on;

  end

end

%% best fit for 6.4
pars3.beta = [ 478.2219  285.2001
  454.2036  585.3093];


pars3.phi =  1.0e-07 *[0.0522    0.0314
    0.1909    0.0192];

pars3.r = [0.4584
    0.4725];

pars3.tau = [5.6604    2.2445
    1.4071    7.3520];

pars3.eta = [ 0.1767    0.4455
    0.7107    0.1360];

pars3.NE = [  83    70
    87    70];

%% best fit whole data

pars3.beta = [ 112.5535  654.1976
  553.3858   89.1597];

pars3.phi = 1.0e-07 *[ 0.0220    0.2074
    0.3414    0.0189];

pars3.r = [0.4457
    0.4576];

pars3.tau = [8.1584    5.5735
    1.8731    8.4054];

pars3.eta = [0.1226    0.1794
    0.5339    0.1190];

pars3.NE = [   83    70
    87    70];


%%


load("./../community/results/psa_seiv_2.mat");
tvec = 0:0.1:15.75;
color_ofthe_fit = [0.8 0.8 0.8];

linewidth = 6;

%%
subplot(2,5,4)
errorbar(time/60,mean(1e3*host4'),std(1e3*host4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue4,Color = blue4);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA H100','FontSize',18);



subplot(2,5,5)
errorbar(time/60,mean(1e3*host5'),std(1e3*host5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue5,Color = blue5);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA 13-15','FontSize',18);



subplot(2,5,9)
errorbar(time/60,mean(1e3*virus4'),std(1e3*virus4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red4,Color = red4);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('PSA HP1','FontSize',18);



subplot(2,5,10)
errorbar(time/60,mean(1e3*virus5'),std(1e3*virus5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red5,Color = red5);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('PSA HS6','FontSize',18);

%%
% pars_median = update_pars(pars2,median(chain_final(5000:10000,:)),mcmcpars2);
% pars_median.eta(pars_median.tau>0) = 1./pars_median.tau(pars_median.tau>0);
% [t2,S_median,V_median,D_median,I_median,E_median] =  simulate_ode(model,pars_median,tvec,pars2.S0,pars2.V0); % mcmc parameter set
% subplot(2,5,4)
% plot(t2,S_median(:,1),'k--','LineWidth',2)
% 
% subplot(2,5,5)
% plot(t2,S_median(:,2),'k--','LineWidth',2)
% 
% subplot(2,5,9)
% plot(t2,V_median(:,1),'k--','LineWidth',2)
% 
% subplot(2,5,10)
% plot(t2,V_median(:,2),'k--','LineWidth',2)

%%
