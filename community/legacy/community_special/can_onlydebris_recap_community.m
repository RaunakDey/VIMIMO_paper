clc;
clear all;

addpath(genpath('./../'))
load('./../revision/combined_posteriors.mat');

load('./../data/triplicate_data.mat');
linewidth = 2;

color_ofthe_fit = [1 0 0]*0.5;
color_ofthe_fit = [0 0 0];
color_ofthe_fill = [0.95 0 0];
transparency = 0.25;
tvec = 0:0.01:15.75; % for better viz
%% new parameters

pars2 = update_pars(pars1,pars_from_dist(chain_stored4(5001:end,:)),mcmcpars);
pars2.epsilon = ones(1,10);
pars2.prob = [0 0 0 0 0]';


%put in dummy value of beta23 as 100

pars2.beta = [0  1.18         0         0         0;
  118.2  60.5  100         0         0;
         0         0   8.2         0         0;
         0         0         0  99  437;
         0         0         0  93  413];


% put in dummy value of tau23 as 1
pars2.tau =[0 1.89 0 0 0;
            1.749 2.19 1 0 0;
            0 0 1.9 0 0 ;
            0 0 0 1.47 2.14;
            0 0 0 1.42 1.97];

pars2.phi = 1e-8 *[         0    1.45        0         0      0;
                             5.14    9.57    5       0         0;
                                  0         0    1.227        0         0;
                                  0         0         0    1.56   6.46;
                                  0         0         0    1.31  8.02];


pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);
pars2.eta

% r
pars2.r = [0.19,0.245,0.22,0.28,0.25]';



%beta_onestep_error = [68.8,2.02,50,0,21.9,43.8,50,79.5,86];
%tau_onestep_error = [0.18,0.56,0.19,0,0.4,0.13,0.12,0.13,0.09];
%phi_onestep_error = [9.57e-9, 1.92e-8,1.39e-8,0,1.68e-8, 2.29e-8,1.99e-8,1.75e-8,1.7e-8];



pars2.Dc = 50e6;
pars2.Dc2 = 30e5;
pars2.Dc3 = 60e5;
pars2.Dc4 = 25e4;
pars2.Dc5 = 15e4;



% NE


pars2.NE = 100*pars2.M;





max_NE = round(max(max(pars2.NE)));
model = SEIVD_diff_NE_diff_debris_abs(5,5,max_NE);
model.host_growth = 0;
model.viral_decay = 0;
model.viral_adsorb = 0;
model.lysis_reset = 0;
model.debris_inhib = 2;
model.debris_inhib2 = 2;
model.debris_inhib3 = 2;
model.debris_inhib4 = 2;
model.debris_inhib5 = 2;

model.diff_beta = 0;


model.name = 'SEIVD-diffabs';


%% going to only take hosts and viruses till a cutoff time.


[t2,S_median,V_median,D_median,I_median,E_median] =  simulate_ode(model,pars2,tvec,pars2.S0,pars2.V0); % mcmc parameter set

%%

hf4 = figure;
subplot(2,5,1)
errorbar(time/60,mean(1e3*host1'),std(1e3*host1'),'o','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255], 'Color',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times


    xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
ylabel({'Host density';'(cells/ml)'});
title('CBA 4','FontSize',18);
    



subplot(2,5,2)
errorbar(time/60,mean(1e3*host2'),std(1e3*host2'),'o','MarkerSize',4,  'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],'Color',[70/255,130/255,180/255] );hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times

   
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 18','FontSize',18);
   




subplot(2,5,3)
errorbar(time/60,mean(1e3*host3'),std(1e3*host3'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],'Color',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times

   
    xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
title('CBA 38','FontSize',18);




subplot(2,5,4)
errorbar(time/60,mean(1e3*host4'),std(1e3*host4'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255], 'Color',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times

  
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA H100','FontSize',18);
    
 




subplot(2,5,5)
errorbar(time/60,mean(1e3*host5'),std(1e3*host5'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],'Color',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times

   
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA 13-15','FontSize',18);
    


%xlabel("Time (hours)");
%ylabel("Host density (cell/ml)");


%%

subplot(2,5,6)
errorbar(time/60,mean(1e3*virus1'),std(1e3*virus1'),'o','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],'Color',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times

  
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e4 1e6 1e8 1e10]);
ylabel({'Phage density';'(virions/ml)'});
title('\phi18:2','FontSize',18);


subplot(2,5,7)
errorbar(time/60,mean(1e3*virus2'),std(1e3*virus2'),'o','MarkerSize',4,  'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],'Color',[70/255,130/255,180/255] );hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times


  xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e4 1e6 1e8 1e10]);
    title('\phi18:3','FontSize',18);
    





subplot(2,5,8)
errorbar(time/60,mean(1e3*virus3'),std(1e3*virus3'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],'Color',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times

  xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('\phi38:1','FontSize',18);


subplot(2,5,9)
errorbar(time/60,mean(1e3*virus4'),std(1e3*virus4'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],'Color',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times

   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('PSA HP1','FontSize',18);




subplot(2,5,10)
errorbar(time/60,mean(1e3*virus5'),std(1e3*virus5'),'o','MarkerSize',4, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],'Color',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times

 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('PSA HS6','FontSize',18);

    %legend('Data','95% confidence interval','Bayesian fit');
    %legend('Box','off');

han=axes(hf4,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
xlabel("Time (hours)");



%% figures




for counts = 1:10000

dummy_beta = 100 + 50*randn;
dummy_tau = 1 + 0.5*randn;



pars2.Dc = 50e6 + 5e6*randn;
pars2.Dc2 = 30e5 + 5e5*randn;
pars2.Dc3 = 60e5 + 5e5*randn;
pars2.Dc4 = 55e4 + 5e4*randn;
pars2.Dc5 = 45e4 + 5e4*randn;



% rewrite pars1
pars2.beta = [0  1.18         0         0         0;
  118.2  60.5  dummy_beta         0         0;
         0         0   8.2         0         0;
         0         0         0  99  437;
         0         0         0  93  413];


% put in dummy value of tau23 as 1
pars2.tau =[0 1.89 0 0 0;
            1.749 2.19 dummy_tau 0 0;
            0 0 1.9 0 0 ;
            0 0 0 1.47 2.14;
            0 0 0 1.42 1.97];

pars2.eta(pars2.tau>0) = 1./pars2.tau(pars2.tau>0);





[t2,S2,V2,~] = simulate_ode(model,pars2,tvec,pars2.S0,pars2.V0); % initial parameter set

%color_of_the_lines = [0.75 0.75 0.75];
%color_of_the_lines = [70/255,210/255,130/255];
%color_of_the_lines = [111,193,157]./255;
%color_of_the_lines=[171,193,157]./255;
color_of_the_lines=[221,93,57]./255;


for i = 1:5
subplot(2,5,i)
patchline(t2,S2(:,i),'edgecolor',color_of_the_lines,'linewidth',2,'edgealpha',0.1);hold on;
ylim([1e4 1e8])
subplot(2,5,i+5)
patchline(t2,V2(:,i),'edgecolor',color_of_the_lines,'linewidth',2,'edgealpha',0.1);hold on;
ylim([1e4 1e11]);

end


end




