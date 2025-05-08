load("./../community/results/psa_seiv_2.mat");
tvec = 0:0.1:15.75;


pars2.NE = [138   156
    108   187];
 model = SEIV_diff_NE(2,2,max(pars2.NE(:)));


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

theta_lhs = [ 342.4616  652.1862  575.1192  489.7728   -7.4812   -8.0501   -7.6421   -8.9611    1.1812    4.1561    2.9987    5.5968    0.5268    0.5722];

theta_optimized = search_minimum_psa(theta_lhs,theta_sd,data,model,pars2,mcmcpars2,0.43,10,1);

theta_optimized_paper = [342.4616  752.1862  475.1192  389.7728   -7.4812   -8.0501   -7.6421   -8.9611    1.3812    2.5561    3.7987    7.5968    0.5268    0.5722];

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



