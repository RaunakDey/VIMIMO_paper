clc;    
clear all;

load('./../community/results/v23.mat');
load('./../community/data/triplicate_data.mat');
%color_of_the_fit = [250 218 221]/255;
color_of_the_fit = [171,193,157]./255;
blue_color_old = [70/255,130/255,180/255];

%%
transient_id = 30000;

close all;

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

transparency = 0.2;
color = [70/255,130/255,180/255];

figure(1)

%%
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

  % that 95% confidence interval limit.
  % you can run for all traces and find the quantile at which that
  % "Half-Gaussian" is cutoff (half-gaussian is an approximation).

  cutoff = 150;
   

count=1;

%%

for i = 30000:100:50000
    i
    
    error_store(count) = loglikefun(chain(i,:),data,pars2,mcmcpars,model,0);
    count=count+1;

    [~, err_vec,~,~,~] = loglikefun(chain(i,:),data,pars2,mcmcpars,model,0);


 pars_samples = update_pars(pars2,chain(i,:),mcmcpars);
 pars_samples.beta2 = pars_samples.beta;

[t3,S3,V3,~] = simulate_ode(model,pars_samples,tvec,pars2.S0,pars2.V0);

    figure(1)

    subplot(2,5,6)
    if err_vec(7) < cutoff
    patchline(t3,V3(:,1),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;
    end

    if err_vec(8) < cutoff
    subplot(2,5,7)
    patchline(t3,V3(:,2),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;
    end


    if err_vec(9) < 25
    subplot(2,5,8)
    patchline(t3,V3(:,3),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;
    end

    if err_vec(10) < 10
    subplot(2,5,9)
    patchline(t3,V3(:,4),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;
    end


    if err_vec(11) < 15
    subplot(2,5,10)
    patchline(t3,V3(:,5),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;
    end


    if err_vec(2) < 10
    subplot(2,5,1)
    patchline(t3,S3(:,1),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;
    end


    if err_vec(3) < 10
    subplot(2,5,2)
    patchline(t3,S3(:,2),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;
    end


    if err_vec(4) < 7
    subplot(2,5,3)
    patchline(t3,S3(:,3),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;
    end


    if err_vec(5) < 7
    subplot(2,5,4)
    patchline(t3,S3(:,4),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;
    end


    if err_vec(6) < 7 
    subplot(2,5,5)
    patchline(t3,S3(:,5),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;
    end


end



%% touch up

figure(1)

subplot(2,5,1)
hold on;
set(gca,'FontSize',18);
title('CBA 4');
ylabel({'Host density','(cells/ml)'})
xticks(0:2:16);
ylim([1e4 1e9]);
yticks([1e4 1e5, 1e6, 1e7, 1e8 1e9]);
xlim([0 16]);
axis square;

subplot(2,5,2)
hold on;
set(gca,'FontSize',18);
title('CBA 18');
xticks(0:2:16);
ylim([1e4 1e9]);
yticks([1e4 1e5, 1e6, 1e7, 1e8 1e9]);
xlim([0 16]);
axis square;

subplot(2,5,3)
hold on;
set(gca,'FontSize',18);
title('CBA 38');
xticks(0:2:16);
ylim([1e4 1e9]);
yticks([1e4 1e5, 1e6, 1e7, 1e8 1e9]);
xlim([0 16]);
axis square;


subplot(2,5,4)
hold on;
set(gca,'FontSize',18);
title('PSA H100');
xticks(0:2:16);
ylim([1e4 1e9]);
yticks([1e4 1e5, 1e6, 1e7, 1e8 1e9]);
xlim([0 16]);
axis square;

subplot(2,5,5)
hold on;
set(gca,'FontSize',18);
title('PSA 13-15');
xticks(0:2:16);
ylim([1e4 1e9]);
yticks([1e4 1e5, 1e6, 1e7, 1e8 1e9]);
xlim([0 16]);
axis square;


subplot(2,5,6)
hold on;
set(gca,'FontSize',18);
title('\phi18:2');
ylabel({'Phage density','(virions/ml)'})
xticks(0:2:16);
ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]);
xlim([0 16]);
axis square;



subplot(2,5,7)
hold on;
set(gca,'FontSize',18);
title('\phi18:3');
xticks(0:2:16);
ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]);
xlim([0 16]);
axis square;


subplot(2,5,8)
hold on;
set(gca,'FontSize',18);
title('\phi38:1');
xticks(0:2:16);
ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]);
xlabel('Time (hrs)');
xlim([0 16]);
axis square;


subplot(2,5,9)
hold on;
set(gca,'FontSize',18);
title('PSA HP1');
xticks(0:2:16);
ylim([1e4 1e11]);
xlim([0 16]);
yticks([1e4, 1e6, 1e8, 1e10]);
axis square;


subplot(2,5,10)
hold on;
set(gca,'FontSize',18);
title('PSA HS6');
xticks(0:2:16);
ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]);
xlim([0 16]);
axis square;







%% Settings for running the model

S0(1) = mean(1e3*host1(1,:)');
S0(2) = mean(1e3*host2(1,:)');
S0(3) = mean(1e3*host3(1,:)');
S0(4) = mean(1e3*host4(1,:)');
S0(5) = mean(1e3*host5(1,:)');


V0(1) = mean(1e3*virus1(1,:)');
V0(2) = mean(1e3*virus2(1,:)');
V0(3) = mean(1e3*virus3(1,:)');
V0(4) = mean(1e3*virus4(1,:)');
V0(5) = mean(1e3*virus5(1,:)');

model =  SEIVD_diff_NE_diff_debris_abs(5,5,200);
model.name = 'SEIVD-diffabs';
model.debris_inhib = 2;
model.debris_inhib2 = 2;
model.debris_inhib3 = 2;
model.debris_inhib4 = 2;
model.debris_inhib5 = 2;


% the important parameters.
pars_afterinf.NH = 5;
pars_afterinf.NV = 5;
pars_afterinf.r = [0.17692
      0.22069
      0.29393
      0.66577
      0.52807];

pars_afterinf.M =[0   1   0   0   0
   1   1   1   0   0
   0   0   1   0   0
   0   0   0   1   1
   0   0   0   1   1];
pars_afterinf.NE = 200*pars_afterinf.M;

pars_afterinf.beta= [         0       2.8875            0            0            0
        194.9       204.58       100.46            0            0
            0            0       19.922            0            0
            0            0            0       525.38       60.653
            0            0            0       488.06       51.326];

pars_afterinf.tau = [  0       2.9978            0            0            0
       1.7394        2.746       2.3189            0            0
            0            0        1.988            0            0
            0            0            0        1.822       4.7139
            0            0            0       2.3157       1.9868];

% BIG PROBLEM HERE.
pars_afterinf.eta= [         0      0.33333            0            0            0
      0.58824      0.37037      0.43478            0            0
            0            0          0.5            0            0
            0            0            0      0.55556      0.21277
            0            0            0      0.43478          0.5];


pars_afterinf.phi = [            0   5.8924e-08            0            0            0
   1.5318e-08   7.7883e-08   2.4067e-08            0            0
            0            0   7.9035e-08            0            0
            0            0            0    6.114e-08   1.1905e-08
            0            0            0   6.0123e-08   2.2635e-08];

%pars_afterinf.phi = zeros(5,5);

pars_afterinf.Dc = 5.0415e+06;
pars_afterinf.Dc2 = 5.9627e+06;
pars_afterinf.Dc3 =  1.1927e+07;
pars_afterinf.Dc4 = 1.8324e+06;
pars_afterinf.Dc5 =  1.52e+06;


pars_afterinf.V0 = [4.2887e+05
   2.8689e+05
     5.28e+05
   1.1033e+05
    1.151e+07];

pars_afterinf.S0 = [2.5111e+06
   5.6423e+06
   3.0257e+06
    6.205e+06
   7.7533e+06];


pars_afterinf.V0 =  V0;
pars_afterinf.S0 =  S0;

pars_afterinf.epsilon = [1 1 1 1 1 1 1 1 1 1];

% not important -- the model class does not call these paramters, the
% plots should not change if these parameters are changed.
pars_afterinf.a = eye(5);
pars_afterinf.m = [  0.00081472
   0.00090579
   0.00012699
   0.00091338
   0.00063236];


pars_afterinf.q = zeros(5,1);
pars_afterinf.prob = [    0
     0
     0
     0
     0];


%%
pars_afterinf.beta= [         0       2.8911            0            0            0
        211.2647      203.5891      98.5805            0            0
            0            0       20.2873            0            0
            0            0            0       525.5791      64.1598
            0            0            0       490.9018    54.9269];

pars_afterinf.phi = 10.^[            0   -7.3932            0            0            0
   -7.9722   -7.2541   -7.8884            0            0
            0            0    -7.1031           0            0
            0            0            0    -7.5251   -7.4906
            0            0            0    -7.2127    -7.7158];
% tbd
pars_afterinf.tau = [  0       2.9091            0            0            0
       1.7654         2.7550    2.4327           0            0
            0            0        2.1160            0            0
            0            0            0        1.7632       4.6215
            0            0            0        2.3340      2.1620];

pars_afterinf.eta = 1./pars_afterinf.tau;

pars_afterinf.eta (isinf(pars_afterinf.eta )) = 0 ;


pars_afterinf.r = [0.1738 
      0.1966
      0.2807
      0.6626
      0.5329];



pars_afterinf.Dc = 0.4908e7;
pars_afterinf.Dc2 = 0.6055e7;
pars_afterinf.Dc3 =  1.0982e7;
pars_afterinf.Dc4 = 0.2175e7;
pars_afterinf.Dc5 =  0.1567e7;

%%
tvec = 0:0.1:15.75;
[t_after,S_after,V_after,D_after,I_after,E_after] =  simulate_ode(model,pars_afterinf,tvec,pars_afterinf.S0,pars_afterinf.V0); % mcmc parameter set



clear ans;
%color_ofthe_fit = [171,230,157]./255;
color_ofthe_fit = [0,0,255]./255;
linewidth = 3;

figure(1);

subplot(2,5,1)
plot(t_after,S_after(:,1),'-','Color',color_ofthe_fit,'LineWidth',linewidth,'LineStyle','--');


subplot(2,5,2)
plot(t_after,S_after(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth,'LineStyle','--');

subplot(2,5,3)
plot(t_after,S_after(:,3),'-','Color',color_ofthe_fit,'LineWidth',linewidth,'LineStyle','--');


subplot(2,5,4)
plot(t_after,S_after(:,4),'-','Color',color_ofthe_fit,'LineWidth',linewidth,'LineStyle','--');


subplot(2,5,5)
plot(t_after,S_after(:,5),'-','Color',color_ofthe_fit,'LineWidth',linewidth,'LineStyle','--');


subplot(2,5,6)
plot(t_after,V_after(:,1),'-','Color',color_ofthe_fit,'LineWidth',linewidth,'LineStyle','--');


subplot(2,5,7)
plot(t_after,V_after(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth,'LineStyle','--');

subplot(2,5,8)
plot(t_after,V_after(:,3),'-','Color',color_ofthe_fit,'LineWidth',linewidth,'LineStyle','--');


subplot(2,5,9)
plot(t_after,V_after(:,4),'-','Color',color_ofthe_fit,'LineWidth',linewidth,'LineStyle','--');


subplot(2,5,10)
plot(t_after,V_after(:,5),'-','Color',color_ofthe_fit,'LineWidth',linewidth,'LineStyle','--');


