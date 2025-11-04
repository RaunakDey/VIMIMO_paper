clear all;
clc;

%% 


color_green = [171,193,157]./255;



blue1 = [179,205,224]/255;
blue2 = [100,151,177]/255;
blue3 = [0,91,150]/255;
blue4 = [3,10,208]/255;
blue5 = [0,0,75]/255;



    

%% community (5H)
load('./../../community/data/nophage_control.mat');

upto = 8;
CBA4 = 1e3*CBA4(1:upto,:);
CBA18 = 1e3*CBA18(1:upto,:);
CBA38 = 1e3*CBA38(1:upto,:);
H100 = 1e3*H100(1:upto,:); 
PSA13_15 = 1e3*PSA13_15(1:upto,:);
time = time(1:upto);

 [slope_ln_cba4,std1]  = get_growth_rate_exponential(time, CBA4);
 [slope_ln_cba18,std2] = get_growth_rate_exponential(time, CBA18);
 [slope_ln_cba38,std3] = get_growth_rate_exponential(time, CBA38);
[slope_ln_h100,std4] = get_growth_rate_exponential(time, H100);
[slope_ln_1315,std5]  = get_growth_rate_exponential(time, PSA13_15);


gr_community = [slope_ln_cba4, slope_ln_cba18, slope_ln_cba38, slope_ln_h100, slope_ln_1315];
std_community = [std1,std2,std3,std4,std5];
%% community (5H 5V)
load('./../../community/data/triplicate_data.mat');

cba4_infected = 1e3*host1(1:upto,:);
cba18_infected = 1e3*host2(1:upto,:);
cba38_infected = 1e3*host3(1:upto,:);
h100_infected = 1e3*host4(1:upto,:);
psa1315_infected = 1e3*host5(1:upto,:);
time = time(1:upto)/60;

 [slope_ln_cba4_inf,std1_inf]  = get_growth_rate_exponential(time, cba4_infected);
 [slope_ln_cba18_inf,std2_inf]  = get_growth_rate_exponential(time, cba18_infected);
 [slope_ln_cba38_inf,std3_inf]  = get_growth_rate_exponential(time, cba38_infected);
[slope_ln_h100_inf,std4_inf]  = get_growth_rate_exponential(time, h100_infected);
[slope_ln_1315_inf,std5_inf]   = get_growth_rate_exponential(time, psa1315_infected);

gr_inf_community = [slope_ln_cba4_inf, slope_ln_cba18_inf, slope_ln_cba38_inf, slope_ln_h100_inf, slope_ln_1315_inf];
std_inf_community = [std1_inf,std2_inf,std3_inf,std4_inf,std5_inf];


%%
gap = 0.1;

%% plots

load('./../../community/data/triplicate_data.mat');
host1 =host1*1e3;
host2 =host2*1e3;
host3 =host3*1e3;
host4 =host4*1e3;
host5 =host5*1e3;
load('./../../community/data/nophage_control.mat');
CBA4 = CBA4*1e3;
CBA18 = CBA18*1e3;
CBA38 = CBA38*1e3;
H100 = H100*1e3;
PSA13_15 = PSA13_15*1e3;

market_size_given = 10;


figure(2)
subplot(2,3,1)
errorbar(time,mean(host1,2),std(host1'),'LineWidth',1,'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',blue1,'MarkeredgeColor','k',Color=blue1 );hold on;
hold on;
errorbar(time,mean(CBA4,2),std(CBA4'), 'LineWidth',1,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k',Color=blue1 )
ylim([1e5 1e9]);
yticks([1e5,1e6, 1e7,1e8,1e9]);
set(gca,'FontSize',18);
set(gca,'YScale','log');
ylabel({'Bacterial density';'(cells/mL)'});
xlabel('Time (hrs)');
xticks(0:2:16);
title('CBA 4');
set(gca,'fontname','times');
xlim([0 16]);
xticks(0:2:16);


subplot(2,3,2)
errorbar(time,mean(host2,2),std(host2'),'LineWidth',1,'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',blue2,'MarkerEdgeColor','k',Color=blue2 );hold on;
hold on;
errorbar(time,mean(CBA18,2),std(CBA18'), 'LineWidth',1,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k',Color=blue2 )
ylim([1e5 1e9]);
yticks([1e5,1e6, 1e7,1e8,1e9]);
set(gca,'FontSize',18);
set(gca,'YScale','log');
xlabel('Time (hrs)');
xticks(0:2:16);
title('CBA 18');
set(gca,'fontname','times');
xlim([0 16]);
xticks(0:2:16);



subplot(2,3,3)
errorbar(time,mean(host3,2),std(host3'),'LineWidth',1,'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',blue3,'MarkerEdgeColor','k',Color=blue3 );hold on;
hold on;
errorbar(time,mean(CBA38,2),std(CBA38'), 'LineWidth',1,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k',Color=blue3 )
ylim([1e5 1e9]);
yticks([1e5,1e6, 1e7,1e8,1e9]);
set(gca,'FontSize',18);
set(gca,'YScale','log');
xlabel('Time (hrs)');
xticks(0:2:16);
title('CBA 38');
set(gca,'fontname','times');
xlim([0 16]);
xticks(0:2:16);



subplot(2,3,4)
errorbar(time,mean(host4,2),std(host4'),'LineWidth',1,'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',blue4,'MarkerEdgeColor','k',Color=blue4 );hold on;
hold on;
errorbar(time,mean(H100,2),std(H100'), 'LineWidth',1,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k',Color=blue4)
ylim([1e5 1e9]);
yticks([1e5,1e6, 1e7,1e8,1e9]);
set(gca,'FontSize',18);
set(gca,'YScale','log');
ylabel({'Bacterial density';'(cells/mL)'});
xlabel('Time (hrs)');
xticks(0:2:16);
title('PSA H100');
set(gca,'fontname','times');
xlim([0 16]);
xticks(0:2:16);



subplot(2,3,5)
errorbar(time,mean(host5,2),std(host5'),'LineWidth',1,'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',blue5,'MarkerEdgeColor','k',Color=blue5);hold on;
hold on;
errorbar(time,mean(PSA13_15,2),std(PSA13_15'), 'LineWidth',1,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k',Color=blue5)
ylim([1e5 1e9]);
yticks([1e5,1e6, 1e7,1e8,1e9]);set(gca,'FontSize',18);
set(gca,'YScale','log');
xlabel('Time (hrs)');
xticks(0:2:16);
title('PSA 13-15');
set(gca,'fontname','times');
xlim([0 16]);
xticks(0:2:16);




subplot(2,3,6)
hold on;
errorbar((1:5)+gap,gr_community,std_community,'LineWidth',1,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k')
errorbar((1:5)-gap,gr_inf_community,std_inf_community,'o','MarkerSize',market_size_given,'MarkerFaceColor','b','Color','b','LineWidth',1);
set(gca,'FontSize',20,'FontName','times');
xticks(1:5);
xticklabels({'CBA 4','CBA 18','CBA 38','PSA H100','PSA 13-15'});
set(gca,'XTickLabelRotation',90)
ylabel({'Growth rates','r (cells/hr)'});
xline(3.5,'--k',LineWidth=2.5);
title({"Exponential fits to the";"first 4 hours of the community"})

legend("Phage-free control","Phage-infected community");
box on;

