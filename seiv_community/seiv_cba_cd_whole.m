clc;
clear all;

load("./../community/results/cba_seiv_2.mat");

tvec = 0:0.1:15.75;



 pars2.NE = [0 10 0
     135   172   100
     0     0    77];
 model = SEIV_diff_NE(3,3,max(pars2.NE(:)));


theta_sd = [10
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

theta_lhs = [368.8264  113.6682  124.1114   63.2691   81.8345   -8.7005   -8.0013   -7.5070   -8.9248   -7.3836    1.9432    3.4407    3.6011    2.7950  3.9504    0.2265    0.1223    0.1594];


%theta_optimized = search_minimum_cba(theta_lhs,theta_sd,data,model,pars2,mcmcpars2,1,15,1)

theta_optimized_paper = [478.8264  203.6682  124.1114  183.2691  191.8345   -8.9005   -8.2013   -7.9070   -8.5248   -7.3836    1.5432    2.5407    3.2011    4.2950    3.8504    0.0765    0.0723    0.0594];
theta_optimized =theta_optimized_paper;

pars2 = update_pars(pars2,theta_optimized,mcmcpars2);

pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);
[t2,S_median,V_median,D_median,I_median,E_median] =  simulate_ode(model,pars2,tvec,pars2.S0,pars2.V0); % mcmc parameter set



%%
hf4 = figure;


subplot(2,3,1)
errorbar(time/60,mean(1e3*host1'),std(1e3*host1'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA H100','FontSize',18);
    
    plot(t2,S_median(:,1),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,3,2)
errorbar(time/60,mean(1e3*host2'),std(1e3*host2'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA 13-15','FontSize',18);
    
    plot(t2,S_median(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth);


    subplot(2,3,3)
errorbar(time/60,mean(1e3*host3'),std(1e3*host3'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA 13-15','FontSize',18);
    
    plot(t2,S_median(:,3),'-','Color',color_ofthe_fit,'LineWidth',linewidth);



subplot(2,3,4)
errorbar(time/60,mean(1e3*virus1'),std(1e3*virus1'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('PSA HP1','FontSize',18);
    
    plot(t2,V_median(:,1),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,3,5)
errorbar(time/60,mean(1e3*virus2'),std(1e3*virus2'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('PSA HS6','FontSize',18);
    
    plot(t2,V_median(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
    

subplot(2,3,6)
errorbar(time/60,mean(1e3*virus3'),std(1e3*virus3'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('PSA HS6','FontSize',18);
    
    plot(t2,V_median(:,3),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
 


han=axes(hf4,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
xlabel("Time (hours)");

