
load('./../community/results/cba_seiv_2.mat');
%%
color_ofthe_fit = [0.8 0.8 0.8];
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

upper_limit = [1000 1000 1000 1000 1000 -7 -7 -7 -7 -7 1 10 10 10 100 0.3 0.3 0.3];
lower_limit = [10 10 10 10 10 -9 -9 -9 -9 -9 0.5 0.5 0.5 0.5 0.5 0.1 0.1 0.1];
draws = 1000;

sampled_params = lhs_sampling(lower_limit,upper_limit, draws);
disp(sampled_params);

%%

transparency = 0.1;
linewidth = 6;
skip = 1;

for i = 1:draws

    i
pars2 = update_pars(pars2,sampled_params(i,:),mcmcpars2);
pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);
[t2,S_median,V_median,D_median,I_median,E_median] =  simulate_ode(model,pars2,tvec,pars2.S0,pars2.V0); % mcmc parameter set

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


%%


upper_limit = [1000 1000 1000 1000 1000 -7 -7 -7 -7 -7 1 10 10 10 100 0.3 0.3 0.3];
lower_limit = [10 10 10 10 10 -9 -9 -9 -9 -9 0.5 0.5 0.5 0.5 0.5 0.1 0.1 0.1];
draws = 1000;

sampled_params = lhs_sampling(lower_limit,upper_limit, draws);
disp(sampled_params);




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
transparency = 0.05;

tvec = 0:0.1:15.75;
for i = 30000:skip:50000

    i
  if error_psa_short(chain_final(i,:),pars2,data,model,mcmcpars2,0.5) ~= Inf
pars2 = update_pars(pars2,chain_final(i,:),mcmcpars2);
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

%%


load("./../community/results/psa_seiv_2.mat");
tvec = 0:0.1:15.75;
color_ofthe_fit = [0.8 0.8 0.8];
linewidth = 6;

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
pars_median = update_pars(pars2,median(chain_final(5000:10000,:)),mcmcpars2);
pars_median.eta(pars_median.tau>0) = 1./pars_median.tau(pars_median.tau>0);
[t2,S_median,V_median,D_median,I_median,E_median] =  simulate_ode(model,pars_median,tvec,pars2.S0,pars2.V0); % mcmc parameter set
subplot(2,5,4)
plot(t2,S_median(:,1),'k--','LineWidth',2)

subplot(2,5,5)
plot(t2,S_median(:,2),'k--','LineWidth',2)

subplot(2,5,9)
plot(t2,V_median(:,1),'k--','LineWidth',2)

subplot(2,5,10)
plot(t2,V_median(:,2),'k--','LineWidth',2)

%%
