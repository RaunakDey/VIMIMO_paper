clear all;
clc;

%% load data


%color_black = [70/255,130/255,180/255];
color_black = [0,0,0];

color_green =[171,193,157]./255;
color_gray =[0.8,0.8,0.8];


skips = 100;
tvec_long = 0:0.01:3.5;

burn = 5000;

load("./../../one-step-experiments/result_replicates/hp1_all13.mat");
hp1_all_beta = chain(burn:end,4);

load("./../../one-step-experiments/result_replicates/HP1_H10011.mat");
hp1_h100_beta = chain(burn:end,4);

load("./../../one-step-experiments/result_replicates/HP1_13-1511.mat");
hp1_1315_beta = chain(burn:end,4);

figure(1)
% subplot(2,4,4)
% smoothHistogram(hp1_all_beta,10,[0 0 0]);
% hold on;
% smoothHistogram(hp1_h100_beta,10,color_black);
% smoothHistogram(hp1_1315_beta,10,'k');
% 
% set(gca,'FontSize',18);
% xlabel('\beta (virions/cell)');
% ylabel('PDF');
% legend('HP1 on all hosts','HP1 on H100','HP1 on 13-15')
% set(gca,'fontname','times');


beta_psa_comm = [511.9 479.26 60.031 49.744];
beta_psa_comm_error =[ 33.787       45.135         5.77       6.3712];

gap = 0.1;


%%
subplot(2,4,4)
errorbar((1:2)-gap, [mean(hp1_h100_beta),mean(hp1_1315_beta)],[std(hp1_h100_beta),std(hp1_1315_beta)],Marker="o", MarkerEdgeColor=color_black,MarkerFaceColor=color_black,Color=color_black,MarkerSize=10,LineStyle="none")
hold on;
yline(mean(hp1_all_beta),'--');
errorbar((1:2)+gap, beta_psa_comm(1:2),beta_psa_comm_error(1:2),Marker="o",MarkerEdgeColor=color_green,MarkerFaceColor=color_green,Color=color_green,MarkerSize=10,LineStyle='none')
errorbar(3,mean(hp1_all_beta),std(hp1_all_beta),Marker="square",MarkerEdgeColor='k',MarkerFaceColor='w',Color='k',MarkerSize=10,LineStyle='none');
smoothHistogram(hp1_h100_beta,10,color_black);
smoothHistogram(hp1_1315_beta,10,'k');
xticks(1:3);
xticklabels({'H100','PSA 13-15','all hosts'})
xtickangle(90)
xlim([0 4]);
set(gca,'FontSize',18);
ylabel('\beta (virions/cell)');

title('Burst size of HP1');
%legend('HP1 on all hosts','HP1 on H100','HP1 on 13-15')
set(gca,'fontname','times');
ylim([0 600]);

%%

[h1,p1,stat1] = kstest2(hp1_all_beta,hp1_h100_beta,'Tail','smaller','Alpha',0.01)

[h2,p2,stat2] = kstest2(hp1_all_beta,hp1_1315_beta,'Tail','smaller','Alpha',0.01)

%% hs6

burn = 5000;

load("./../../one-step-experiments/result_replicates/hs6_all13.mat");
hs6_all_beta = chain(burn:end,4);

%17
%load("./../../one-step-experiments/result_replicates/hs6_h100_11.mat");
load("./../../one-step-experiments/result_replicates/hs6_h100_18.mat");
% can use 18 if needed to match with the other figure
hs6_h100_beta = chain(burn:end,4);
mean(hs6_h100_beta)

%load("./results_same_phi/HS6_H100_3-inferred.mat");
%hs6_h100_beta = (chain(burn:end,4)) - mean((chain(burn:end,4))) + 391;


%hs6_h100_beta = 0.4*(chain(100:end,4) - mean(chain(100:end,4))) + mean(chain(100:end,4)) ;
%hs6_h100_beta = 0.5*(chain(burn:end,4) - mean(chain(burn:end,4)))+mean(chain(burn:end,4));


load("./../../one-step-experiments/result_replicates/HS6_13-1511.mat");
hs6_1315_beta = chain(burn:end,4);


% subplot(2,4,8)
% smoothHistogram(hs6_all_beta,10,[0 0 0]);
% hold on;
% smoothHistogram(hs6_h100_beta,10,color_black);
% smoothHistogram(hs6_1315_beta,10,'k');
% 
% set(gca,'FontSize',18);
% xlabel('\beta (virions/cell)');
% ylabel('PDF');
% legend('HS6 on all hosts','HS6 on H100','HS6 on 13-15')
% ylim([0 7e-3]);
% set(gca,'fontname','times');
% 

%%

subplot(2,4,8)
errorbar((1:2)-gap, [mean(hs6_h100_beta),mean(hs6_1315_beta)],[std(hs6_h100_beta),std(hs6_1315_beta)],Marker="o",MarkerEdgeColor=color_black,MarkerFaceColor=color_black,Color=color_black,MarkerSize=10,LineStyle="none")
hold on;
errorbar(3,mean(hs6_all_beta),std(hs6_all_beta),Marker="square",MarkerEdgeColor='k',MarkerFaceColor='w',Color='k',MarkerSize=10,LineStyle='none');
errorbar((1:2)+gap, beta_psa_comm(3:4),beta_psa_comm_error(3:4),Marker="o",MarkerEdgeColor=color_green,MarkerFaceColor=color_green,Color=color_green,MarkerSize=10,LineStyle='none')
yline(mean(hs6_all_beta),'--');
xticks(1:3);
xticklabels({'H100','PSA 13-15','all hosts'})
xtickangle(90)
xlim([0 4]);

set(gca,'FontSize',18);
ylabel('\beta (virions/cell)');
title('Burst size of HS6');
%legend('HP1 on all hosts','HP1 on H100','HP1 on 13-15')
set(gca,'fontname','times');
ylim([0 600]);
%%

[h3,p3,stat3] = kstest2(hs6_all_beta,hs6_h100_beta,'Tail','larger','Alpha',0.01)

[h4,p4,stat4] = kstest2(hs6_all_beta,hs6_1315_beta,'Tail','larger','Alpha',0.01)

%% z test

z1 = (mean(hp1_all_beta) - mean(hp1_h100_beta))./sqrt(var(hp1_all_beta) +var(hp1_h100_beta) )
z2 = (mean(hp1_all_beta) - mean(hp1_1315_beta))./sqrt(var(hp1_all_beta) +var(hp1_1315_beta) )

z3 = (mean(hs6_all_beta) - mean(hs6_h100_beta))./sqrt(var(hs6_all_beta) +var(hs6_h100_beta) )
z4 = (mean(hs6_all_beta) - mean(hs6_1315_beta))./sqrt(var(hs6_all_beta) +var(hs6_1315_beta) )



%% plotting



subplot(2,4,1)
%load("./results_same_phi/HP1_H100_6-inferred.mat");
load("./../../one-step-experiments/result_replicates/hp1_all13.mat");
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





%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_black,'LineWidth',2,'MarkerEdgeColor',color_black,'MarkerFaceColor',color_black)

patchline(time2,y_series2(end,:),'edgecolor',color_gray,...
    'linewidth',4,'edgealpha',0.2);hold on;


for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'square','LineStyle','none','Color', color_gray,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','w'); hold on;
    end
end

end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','k','MarkerFaceColor','k','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('PSA HP1—all hosts');
ylim([1e2 1e8]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
%clear max;
xlim([0 2.5]);
%plot(0,S0/100,'^',Color=[0 0 0],MarkerSize=10,MarkerEdgeColor='k',MarkerFaceColor='k')

xlabel("Time (hrs)");
ylabel("Free phage density (Virions/ml)");
set(gca,'fontname','times');


subplot(2,4,2)
%load("./results_same_phi/HP1_H100_6-inferred.mat");
load("./../../one-step-experiments/result_replicates/HP1_H10011.mat");
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
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_black,'LineWidth',2,'MarkerEdgeColor',color_black,'MarkerFaceColor',color_black); hold on;
    end
end


%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_black,'LineWidth',2,'MarkerEdgeColor',color_black,'MarkerFaceColor',color_black)

patchline(time2,y_series2(end,:),'edgecolor',color_gray,...
    'linewidth',4,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','k','MarkerFaceColor','k','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('PSA HP1—PSA H100');
ylim([1e2 1e8]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
%clear max;
xlim([0 2.5]);
%plot(0,S0/100,'^',Color=[0 0 0],MarkerSize=10,MarkerEdgeColor='k',MarkerFaceColor='k')
xlabel("Time (hrs)");
ylabel("Free phage density (Virions/ml)");
set(gca,'fontname','times');


subplot(2,4,3)
%load("./results_same_phi/HP1_13-15_3-inferred.mat");
load("./../../one-step-experiments/result_replicates/HP1_13-1511.mat")
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

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_black,'LineWidth',2,'MarkerEdgeColor',color_black,'MarkerFaceColor',color_black)

for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', 'k','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k'); hold on;
    end
end

patchline(time2,y_series2(end,:),'edgecolor',color_gray,...
    'linewidth',4,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','k','MarkerFaceColor','k','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('PSA HP1—PSA 13-15');
ylim([1e2 1e8]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
%clear max;
xlim([0 1.75]);
%plot(0,S0/100,'^',Color=[0 0 0],MarkerSize=10,MarkerEdgeColor='k',MarkerFaceColor='k')
xlabel("Time (hrs)");
ylabel("Free phage density (Virions/ml)");
set(gca,'fontname','times');





subplot(2,4,5)
%load("./results_same_phi/HP1_H100_6-inferred.mat");
load("./../../one-step-experiments/result_replicates/hs6_all13.mat");
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



%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_black,'LineWidth',2,'MarkerEdgeColor',color_black,'MarkerFaceColor',color_black)

patchline(time2,y_series2(end,:),'edgecolor',color_gray,...
    'linewidth',4,'edgealpha',0.2);hold on;



for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'square','LineStyle','none','Color', color_gray,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','w'); hold on;
    end
end

end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','k','MarkerFaceColor','k','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('PSA HS6—all hosts');
ylim([1e2 1e8]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
%clear max;
xlim([0 2.5]);
%plot(0,S0/100,'^',Color=[0 0 0],MarkerSize=10,MarkerEdgeColor='k',MarkerFaceColor='k')
xlabel("Time (hrs)");
ylabel("Free phage density (Virions/ml)");
set(gca,'fontname','times');



subplot(2,4,6)
%load("./results_same_phi/HS6_H100_3-inferred.mat");

load("./../../one-step-experiments/result_replicates/hs6_h100_15.mat")
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

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_black,'LineWidth',2,'MarkerEdgeColor',color_black,'MarkerFaceColor',color_black)

 for k = 1:num_replicates
     for j = 1:3
     plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', color_black,'LineWidth',2,'MarkerEdgeColor',color_black,'MarkerFaceColor',color_black); hold on;
     end
 end

patchline(time2,y_series2(end,:),'edgecolor',color_gray,...
    'linewidth',4,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','k','MarkerFaceColor','k','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('PSA HS6—PSA H100');
ylim([1e2 1e8]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
%clear max;
xlim([0 3.5]); 
%plot(0,S0/100,'^',Color=[0 0 0],MarkerSize=10,MarkerEdgeColor='k',MarkerFaceColor='k')
xlabel("Time (hrs)");
ylabel("Free phage density (Virions/ml)");
set(gca,'fontname','times');



subplot(2,4,7)
%load("./results_same_phi/HS6_13-15_3-inferred.mat");

load("./../../one-step-experiments/result_replicates/HS6_13-1511.mat")
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

%errorbar(data.xdata, free_phages_mean , std(free_phages')', 'o','LineStyle','none','Color', color_black,'LineWidth',2,'MarkerEdgeColor',color_black,'MarkerFaceColor',color_black)

for k = 1:num_replicates
    for j = 1:3
    plot(data.xdata{k}./60,data.ydata{k}(:,j),'o','LineStyle','none','Color', 'k','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k'); hold on;
    end
end


patchline(time2,y_series2(end,:),'edgecolor',color_gray,...
    'linewidth',4,'edgealpha',0.2);hold on;


end
xtickangle(45)
%plot(0,y_series2(1,1),'MarkerEdgeColor','k','MarkerFaceColor','k','Marker','square','MarkerSize',8);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
title('PSA HS6—PSA 13-15');
ylim([1e2 1e8]);
%xlim([0 3.5]);
xticks([0 0.5 1 1.5 2 2.5 3 3.5]);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box on;
%clear max;
xlim([0 3.0]);
%plot(0,S0/100,'^',Color=[0 0 0],MarkerSize=10,MarkerEdgeColor='k',MarkerFaceColor='k')
xlabel("Time (hrs)");
ylabel("Free phage density (Virions/ml)");
set(gca,'fontname','times');

