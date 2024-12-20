clear all;
clc;


% Please not fig 1a and 1b are schematics not plots.

%% Figure 1C

color_blue = [70/255,130/255,180/255];
color_green =[171,193,157]./255;

% I am changing it to black and grey
color_blue = [0,0,0];
color_green = [0.8,0.8,0.8];

skips = 10;

tvec_long = 0:0.01:3.5;


figure(1)
subplot(2,4,1)
%load("./results_same_phi/CBA18-2_18_5-inferred.mat");
load("./../one-step-experiments/result_replicates/18:2_18_12_13.mat");
chain_effective = chain(5000:end,:);
S0 = mean(S0_replicates);
V0 = mean(V0_replicates);

for i=1:skips:length(chain_effective)

clear y
clear y0
NE =  round(chain_effective(i,5));
y(1,1) = S0;
y(1,2:NE+2) = 0;
y(1,NE+3) = V0;

[time2,y_series2,time_abs,pre_dil] = one_step_simulate(tvec_long,y,chain_effective(i,:),NE, dilution_factor);
%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',6,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)
%plot(time_free_phages,free_phages,'o','LineStyle','none','Color', color_blue,'LineWidth',6,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)

for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',1,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
    end
end


patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'LineWidth',6,'edgealpha',0.2);hold on;

end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18)
title('\phi18:2—CBA 18');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
clear max;
xlim([0 max(data.xdata{1})/60]);
ylim([1e2 1e8])



subplot(2,4,2)
%load("./results_same_phi/CBA18-3_4_10-inferred.mat");
load("./../one-step-experiments/result_replicates/18:3_4_12.mat")
chain_effective = chain(6000:end,:);
S0 = mean(S0_replicates);
V0 = mean(V0_replicates);


for i=1:skips:length(chain_effective)

clear y
clear y0
NE =  round(chain_effective(i,5));
y(1,1) = S0;
y(1,2:NE+2) = 0;
y(1,NE+3) = V0;

[time2,y_series2,time_abs,pre_dil] = one_step_simulate(tvec_long,y,chain_effective(i,:),NE, dilution_factor);

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',6,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)
%plot(time_free_phages,free_phages,'o','LineStyle','none','Color', color_blue,'LineWidth',6,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)



for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',1,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
    end
end


patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'LineWidth',6,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18)
title('\phi18:3—CBA 4');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
clear max;
xlim([0 max(data.xdata{1})/60 ]);
ylim([1e2 1e8])


subplot(2,4,3)
%load("./results_same_phi/CBA18-3_18_7-inferred.mat");
load("./../one-step-experiments/result_replicates/CBA18-3_18_2024_12.mat")
chain_effective = chain(5000:end,:);
S0 = mean(S0_replicates);
V0 = mean(V0_replicates);


for i=1:skips:length(chain_effective)

clear y
clear y0
NE =  round(chain_effective(i,5));
y(1,1) = S0;
y(1,2:NE+2) = 0;
y(1,NE+3) = V0;

[time2,y_series2,time_abs,pre_dil] = one_step_simulate(tvec_long,y,chain_effective(i,:),NE, dilution_factor);

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',6,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)
%plot(time_free_phages,free_phages,'o','LineStyle','none','Color', color_blue,'LineWidth',6,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)


for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',1,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
    end
end

patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'LineWidth',6,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18)
title('\phi18:3—CBA 18');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
clear max;
xlim([0 max(data.xdata{1})/60 ]);
ylim([1e2 1e8])



subplot(2,4,4)
%load("./results_same_phi/CBA38-1_38_6-inferred.mat");
load("./../one-step-experiments/result_replicates/CBA38-1_38_13.mat")
chain_effective = chain(5000:end,:);
S0 = mean(S0_replicates);
V0 = mean(V0_replicates);


for i=1:skips:length(chain_effective)

clear y
clear y0
NE =  round(chain_effective(i,5));
y(1,1) = S0;
y(1,2:NE+2) = 0;
y(1,NE+3) = V0;

[time2,y_series2,time_abs,pre_dil] = one_step_simulate(tvec_long,y,chain_effective(i,:),NE, dilution_factor);

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',6,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)
%plot(time_free_phages,free_phages,'o','LineStyle','none','Color', color_blue,'LineWidth',6,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)

for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',1,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
    end
end

patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'LineWidth',6,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18)
title('\phi38:1—CBA 38');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
clear max;
xlim([0 max(data.xdata{1})/60 ]);
ylim([1e2 1e8])





subplot(2,4,5)
%load("./results_same_phi/HP1_H100_6-inferred.mat");
load("./../one-step-experiments/result_replicates/HP1_H10011.mat");
S0 = mean(S0_replicates);
V0 = mean(V0_replicates);
chain_effective = chain(5000:end,:);


for i=1:skips:length(chain_effective)

clear y
clear y0
NE =  round(chain_effective(i,5));
y(1,1) = S0;
y(1,2:NE+2) = 0;
y(1,NE+3) = V0;

[time2,y_series2,time_abs,pre_dil] = one_step_simulate(tvec_long,y,chain_effective(i,:),NE, dilution_factor);


for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',1,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
    end
end


%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',6,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)

patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'LineWidth',6,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18)
title('PSA HP1—PSA H100');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
%clear max;
xlim([0 2.5]);
ylim([1e2 1e8])





subplot(2,4,6)
%load("./results_same_phi/HP1_13-15_3-inferred.mat");
load("./../one-step-experiments/result_replicates/HP1_13-1511.mat")
S0 = mean(S0_replicates);
V0 = mean(V0_replicates);
chain_effective = chain(5000:end,:);

for i=1:skips:length(chain_effective)

clear y
clear y0
NE =  round(chain_effective(i,5));


y(1,1) = S0;
y(1,2:NE+2) = 0;
y(1,NE+3) = V0;

[time2,y_series2,time_abs,pre_dil] = one_step_simulate(tvec_long,y,chain_effective(i,:),NE, dilution_factor);

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',6,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)

for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',1,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
    end
end

patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'LineWidth',6,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18)
title('PSA HP1—PSA 13-15');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
%clear max;
xlim([0 1.75]);
ylim([1e2 1e8])

%%

subplot(2,4,7)
%load("./results_same_phi/HS6_H100_3-inferred.mat");

%load("./../one-step-experiments/result_replicates/hs6_h100_15.mat")
load("./../one-step-experiments/result_replicates/hs6_h100_18.mat")
% in this run I moved away from the OSU prior -- that was terrible.

S0 = mean(S0_replicates);
V0 = mean(V0_replicates);
chain_effective = chain(5000:end,:);

for i=1:skips:length(chain_effective)

clear y
clear y0
NE =  round(chain_effective(i,5));
y(1,1) = S0;
y(1,2:NE+2) = 0;
y(1,NE+3) = V0;

[time2,y_series2,time_abs,pre_dil] = one_step_simulate(tvec_long,y,chain_effective(i,:),NE, dilution_factor);

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',6,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)

 for k = 1:num_replicates
     for j = 1:3
     plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',1,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
     end
 end

patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'LineWidth',6,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18)
title('PSA HS6—PSA H100');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
%clear max;
xlim([0 3.5]);
ylim([1e2 2e8])   

%%
subplot(2,4,8)
%load("./results_same_phi/HS6_13-15_3-inferred.mat");

load("./../one-step-experiments/result_replicates/HS6_13-1511.mat")

S0 = mean(S0_replicates);
V0 = mean(V0_replicates);
chain_effective = chain(5000:end,:);


for i=1:skips:length(chain_effective)

clear y
clear y0
NE =  round(chain_effective(i,5));
y(1,1) = S0;
y(1,2:NE+2) = 0;
y(1,NE+3) = V0;

[time2,y_series2,time_abs,pre_dil] = one_step_simulate(tvec_long,y,chain_effective(i,:),NE, dilution_factor);

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',6,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)

for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',1,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
    end
end


patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'LineWidth',6,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontName','Times','FontSize',18)
title('PSA HS6—PSA 13-15');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
%clear max;
xlim([0 3.0]);
ylim([1e2 1e8])


%% Figure 1D



x=1:8;


beta_onestep_conventional = [91.2,0.94,27.29,10.5,54.9,54.2,205.8,276.7];
tau_onestep_conventional = [1.25,1,1.41,1,0.75,0.667,1.33,1.25];
phi_onestep_conventional = [1.88e-8,1.83e-7,1.33e-7,1.65e-8,9.9e-8,1.79e-7,1.87e-7,7.6e-8,1e-7];

beta_onestep_conventional_error = [32.91,0.15,6.28,3.18,21.8,15.74,92.7,110.4];
phi_one_step_conventional_error = [3.4e-09     3.99e-08     2.66e-08     4.96e-10     1.37e-08     3.33e-08     3.79e-08     8.64e-09     3.55e-09]';



%beta_onestep = [   126.55        1.901       63.884        36.55       75.092       87.229       431.35       324.15];
%beta_onestep = [   126.55        1.901       63.884        36.55       75.092       87.229      357.3      324.15];
beta_onestep = [   126.55        1.901       63.884        36.55       75.092       87.229      435.6      324.15]; %18 (hs6_h100)


tau_onestep = [   1.8781       1.9068       2.0352       1.8012       1.4589       1.3709       2.22       2.0193];
phi_onestep = [ 1.776e-08   5.6449e-08   1.2386e-07    3.506e-08   1.2251e-07   8.5023e-08    5.6039e-08    4.4429e-08];



beta_onestep_error = [    36.508      0.12939       48.102       45.319       30.193       42.213       100       71.87];
tau_onestep_error = [   0.038254      0.26513      0.17063      0.27281      0.11118     0.075967     0.105     0.078438];
phi_onestep_error = [ 5.2604e-09   1.1261e-08   2.7217e-08   1.2592e-08   1.3211e-08    9.791e-09   7.6718e-09     7.4481e-09];


 r_bayesian = [      0.2531      0.18702      0.23981      0.21843      0.28019      0.25393       0.29      0.25793];
 r_bayesian_std = [  0.013507     0.048525     0.048     0.046714     0.042521     0.049102     0.042063     0.019837];



r_onestep = [r_bayesian(2), mean([r_bayesian(2),r_bayesian(3)]), r_bayesian(4), mean([r_bayesian(5),r_bayesian(7)]),mean([r_bayesian(6),r_bayesian(8)]) ];
r_onestep_error = [r_bayesian_std(2), sqrt(r_bayesian_std(2)^2 + r_bayesian_std(3)^2), r_bayesian_std(4),  sqrt(r_bayesian_std(5)^2 + r_bayesian_std(7)^2), sqrt(r_bayesian_std(6)^2 + r_bayesian_std(8)^2) ];


r_onestep_conventional = [0.19,0.245,0.22,0.28,0.25];
r_onestep_conventional_error =   [0.019093     0.021355     0.030257     0.052883     0.044104]';


boxes_onestep = [   148.01           10       171.79       77.868        138.8       108.48          156     187.61];
boxes_onestep_error = [ 22.555     0.043257       84.682       64.002       72.778       25.955       60       61.193];


color_def2 = [0.9290 0.6940 0.1250];
%color_def = [111,193,157]./255;
color_def = [0,0,0];


figure(2)
subplot(1,5,1)
errorbar(x-0.1,beta_onestep_conventional,beta_onestep_conventional_error,'^','MarkerSize',14,'MarkerFaceColor',color_def2,'Color',color_def2,LineWidth=2);
hold on;
errorbar(x+0.1,beta_onestep,beta_onestep_error,'bo','MarkerSize',14,'MarkerFaceColor',[0.5,0.5,0.5],'Color',[0.5,0.5,0.5],LineWidth=2);
title({'Burst size', '\beta (virions/cell)'});
ylim([1 600]);
set(gca,'FontName','Times','FontSize',18);
xticks(1:8)
xticklabels([1,2,3,5,6,7,8,9]);
yticks(0:100:600);
ylim([-1 600])
xlim([0 9]);


subplot(1,5,2)
plot(x-0.1,tau_onestep_conventional,'^','MarkerSize',14,'MarkerFaceColor',color_def2,'Color',color_def2,LineWidth=2);
hold on;
errorbar(x+0.1,tau_onestep,tau_onestep_error,'bo','MarkerSize',14,'MarkerFaceColor',[0.5,0.5,0.5],'Color',[0.5,0.5,0.5],LineWidth=2);
title({'Latent period',' \tau (hr)'});
ylim([0.2 3]);
set(gca,'FontName','Times','FontSize',18);
xticks(1:8);
xticklabels([1,2,3,5,6,7,8,9]);
xlim([0 9]);


subplot(1,5,4)
errorbar([1 2 3 4 5 6 7 8 9]-0.1,phi_onestep_conventional,phi_one_step_conventional_error,'^','MarkerSize',14,'MarkerFaceColor',color_def2,'Color',color_def2,LineWidth=2);
hold on;
errorbar([1 2 3 5 6 7 8 9]+0.1,phi_onestep,phi_onestep_error,'bo','MarkerSize',14,'MarkerFaceColor',[0.5,0.5,0.5],'Color',[0.5,0.5,0.5],LineWidth=2);
title({'Adsorption rate','\phi (ml/hr)'})
ylim([1e-10 2.5e-7]);
set(gca,'FontName','Times','FontSize',18);
xticks(1:9)
xticklabels([1,2,3,4,5,6,7,8,9]);
xlim([0 10]);
%legend('Conventional','Bayesian');
%legend('boxoff');



cv_mean_onestep  = 1./sqrt(boxes_onestep);
cv_onestep_error = 0.5* (boxes_onestep).^(-1.5).*boxes_onestep_error;
%cv_onestep_error  = 1./sqrt(boxes_onestep + boxes_onestep_error) - 1./sqrt(boxes_onestep + boxes_onestep_error)

subplot(1,5,3)
errorbar(x,cv_mean_onestep,cv_onestep_error,'ko','MarkerSize',14,'MarkerFaceColor',[0.5,0.5,0.5],'Color',[0.5,0.5,0.5],LineWidth=2);
title({'Coefficient of',' variation of \tau (CV)'});
set(gca,'FontName','Times','FontSize',18);
xticks(1:8);
xticklabels([1,2,3,5,6,7,8,9]);
ylim([0 0.35])
%yticks(0:0.1:0.3);
xlim([0 9]);


subplot(1,5,5)
errorbar([1:5] ,r_onestep_conventional,r_onestep_conventional_error,'^','MarkerSize',14,'MarkerFaceColor',color_def2,'Color',color_def2,LineWidth=2);
hold on;
%errorbar([1:5]+0.1, r_onestep, r_onestep_error,'bo','MarkerSize',14,'MarkerFaceColor',[0.5,0.5,0.5],'Color',[0.5,0.5,0.5],LineWidth=2);

title({'Growth rate','r (cells/hr)'});
ylim([0.1 0.35]);
set(gca,'FontName','Times','FontSize',18);
xlim([0.5 5.5])
xticks(1:5);
xticklabels({'CBA 4','CBA 18','CBA 38','PSA H100','PSA 13-15'});
set(gca,'XTickLabelRotation',90);




