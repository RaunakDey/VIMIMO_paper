clear all;
clc;


%load("./results_same_phi/CBA18-2_18_5-inferred.mat");
%load("./results_same_phi/CBA18-3_4_10-inferred.mat");
%load("./results_same_phi/CBA18-3_18_7-inferred.mat");
%load("./results_same_phi/CBA38-1_38_6-inferred.mat");
%load("./results_same_phi/HP1_H100_6-inferred.mat");
%load("./results_same_phi/HP1_13-15_3-inferred.mat");
%load("./results_same_phi/HS6_H100_3-inferred.mat");
%load("./results_same_phi/HS6_13-15_3-inferred.mat");



color_blue = [70/255,130/255,180/255];



color_green =[171,193,157]./255;
skips = 100;

tvec_long = 0:0.01:3.5;



subplot(2,4,1)
%load("./results_same_phi/CBA18-2_18_5-inferred.mat");
load("./result_replicates/18:2_18_12_13.mat");
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
%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)
%plot(time_free_phages,free_phages,'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)

for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
    end
end


patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'linewidth',2,'edgealpha',0.2);hold on;

end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('\phi18:2--CBA 18');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
clear max;
xlim([0 max(data.xdata{1})/60]);


subplot(2,4,2)
%load("./results_same_phi/CBA18-3_4_10-inferred.mat");
load("./result_replicates/18:3_4_12.mat")
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

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)
%plot(time_free_phages,free_phages,'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)



for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
    end
end


patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'linewidth',2,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('\phi18:3--CBA 4');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
clear max;
xlim([0 max(data.xdata{1})/60 ]);


subplot(2,4,3)
%load("./results_same_phi/CBA18-3_18_7-inferred.mat");
load("./result_replicates/CBA18-3_18_2024_12.mat")
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

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)
%plot(time_free_phages,free_phages,'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)


for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
    end
end

patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'linewidth',2,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('\phi18:3--CBA 18');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
clear max;
xlim([0 max(data.xdata{1})/60 ]);
   

subplot(2,4,4)
%load("./results_same_phi/CBA38-1_38_6-inferred.mat");
load("./result_replicates/CBA38-1_38_13.mat")
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

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)
%plot(time_free_phages,free_phages,'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)

for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
    end
end

patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'linewidth',2,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('\phi38:1--CBA 38');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
clear max;
xlim([0 max(data.xdata{1})/60 ]);









subplot(2,4,5)
%load("./results_same_phi/HP1_H100_6-inferred.mat");
load("./result_replicates/HP1_H10011.mat");
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
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
    end
end


%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)

patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'linewidth',2,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('PSA HP1--PSA H100');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
%clear max;
xlim([0 2.5]);





subplot(2,4,6)
%load("./results_same_phi/HP1_13-15_3-inferred.mat");
load("./result_replicates/HP1_13-1511.mat")
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

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)

for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
    end
end

patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'linewidth',2,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('PSA HP1--PSA 13-15');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
%clear max;
xlim([0 1.75]);


subplot(2,4,7)
%load("./results_same_phi/HS6_H100_3-inferred.mat");

load("./result_replicates/hs6_h100_15.mat")
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

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)

 for k = 1:num_replicates
     for j = 1:3
     plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
     end
 end

patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'linewidth',2,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('PSA HS6--PSA H100');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
%clear max;
xlim([0 3.5]);
   


subplot(2,4,8)
%load("./results_same_phi/HS6_13-15_3-inferred.mat");

load("./result_replicates/HS6_13-1511.mat")
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

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)

for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue); hold on;
    end
end


patchline(time2,y_series2(end,:),'edgecolor',color_green,...
    'linewidth',2,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('PSA HS6--PSA 13-15');
%ylim([1e2 1e9]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
%clear max;
xlim([0 3.0]);

%%

% figure(2)
% load('./results_same_phi/HP1_allhosts.mat');
% 
% chain_effective = chain(2000:end,:);
% 
% for i=1:skips:length(chain_effective)
% 
% clear y
% clear y0
% NE =  round(chain_effective(i,5));
% y(1,1) = S0;
% y(1,2:NE+2) = 0;
% y(1,NE+3) = V0;
% 
% [time2,y_series2,time_abs,pre_dil] = one_step_simulate(time_free_phages,y,chain_effective(i,:),NE, dilution_factor);
% 
% errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_blue,'LineWidth',2,'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)
% patchline(time2,y_series2(end,:),'edgecolor',color_green,...
%     'linewidth',2,'edgealpha',0.2);hold on;
% 
% 
% end
% 
% plot(0,y_series2(1,1),'MarkerEdgeColor','r','MarkerFaceColor','r','Marker','square','MarkerSize',8);
% set(gca, 'YScale', 'log');
% set(gca,'FontSize',18)
% title('PSA HP1--all hosts');
% ylim([1e2 1e9]);
% xlim([0 3.5]);
% xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
% yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
% box on;