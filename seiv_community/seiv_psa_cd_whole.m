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


%%

upper_limit = [700 700 700 700 -7 -7 -7 -7 10 10 10 10 0.7 0.7];
lower_limit = [50 50 50 50 -9 -9 -9 -9 0.5 0.5 0.5 0.5 0.3 0.3];
 
draws = 500;
sampled_params = lhs_sampling(lower_limit,upper_limit, draws);
error_psa = [];

for i = 1:draws
    i
    error_psa(i) = error_psa_short(sampled_params(i,:),pars2,data,model,mcmcpars2,1);
    error_psa(i)
end

[value,index] = min(error_psa);
theta_lhs = sampled_params(index,:);

%%

theta_lhs = [177.3423  588.1937  491.9885  501.9945   -7.0920   -7.8204   -7.4877   -8.3562 0.8436    6.1534    9.8564    2.5685    0.4716    0.4648 ];

%theta_optimized = search_minimum_psa(theta_lhs,theta_sd,data,model,pars2,mcmcpars2,1,10,1);

theta_optimized_paper = [  227.3423  888.1937  421.9885  761.9945   -7.0920   -7.8204   -7.4877   -8.3562    1.0436    6.3534   10.6564    2.5685    0.4716    0.4648]; 

pars2 = update_pars(pars2,theta_optimized_paper,mcmcpars2);

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
   title('PSA-HP1','FontSize',18);
    
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
  title('PSA-HS6','FontSize',18);
    
    plot(t2,V_median(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
    %legend('Data','95% confidence interval','Bayesian fit');
    %legend('Box','off');

han=axes(hf4,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
%xlabel("Time (hours)");



