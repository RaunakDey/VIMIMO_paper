clc;
clear all;
addpath(genpath(pwd)); % add current directory to path (mcmcstat is included)

%% no phage control -- 1 strain in a flask


cba4_nophage =  1e3*[2.11E+03	1.30E+03	3.00E+03;
2.71E+03	2.39E+03	2.21E+03;
5.50E+04	7.81E+04	7.20E+04];

cba18_nophage = 1e3*[4.89E+03	7.92E+03	9.53E+03;
1.56E+04	9.52E+03	1.26E+04;
3.71E+05	2.39E+05	1.93E+05];

cba38_nophage  = 1e3*[7.81E+03	1.01E+04	7.77E+03;
1.42E+04	1.04E+04	1.29E+04;
1.01E+05	2.02E+05	2.40E+05];

h100_nophage = 1e3*[1.21E+04	1.24E+04	1.18E+04;
7.36E+04	5.34E+04	6.42E+04;
6.52E+05	6.30E+05	9.31E+05];

psa1315_nophage = 1e3*[1.43E+04	1.56E+04	1.84E+04;
2.86E+04	2.31E+04	2.57E+04;
2.95E+05	2.02E+05	3.21E+05];

time2 = [0,175,945]/60;


%% isolated strains
slope_ln_cba4  = get_growth_rate_exponential(time2, cba4_nophage);
 slope_ln_cba18 = get_growth_rate_exponential(time2, cba18_nophage);
 slope_ln_cba38 = get_growth_rate_exponential(time2, cba38_nophage);
 slope_ln_h100  = get_growth_rate_exponential(time2, h100_nophage);
slope_ln_1315  = get_growth_rate_exponential(time2, psa1315_nophage);

gr_isolated = [slope_ln_cba4, slope_ln_cba18, slope_ln_cba38, slope_ln_h100, slope_ln_1315];

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
market_size_given = 10;
color_green = [171,193,157]./255;
% 
% r_bayesian = [      0.2531      0.18702      0.23981      0.21843      0.28019      0.25393       0.27987      0.25793];
% r_bayesian_std = [  0.013507     0.048525     0.048     0.046714     0.042521     0.049102     0.040756     0.019837];
% % 
% r_onestep = [r_bayesian(2), mean([r_bayesian(2),r_bayesian(3)]), r_bayesian(4), mean([r_bayesian(5),r_bayesian(7)]),mean([r_bayesian(6),r_bayesian(8)]) ];
% r_onestep_error = [r_bayesian_std(2), sqrt(r_bayesian_std(2)^2 + r_bayesian_std(3)^2), r_bayesian_std(4),  sqrt(r_bayesian_std(5)^2 + r_bayesian_std(7)^2), sqrt(r_bayesian_std(6)^2 + r_bayesian_std(8)^2) ];


%%


figure()
hold on;
% uninfected isolated
%plot((1:5)-gap,gr_isolated,'MarkerSize',market_size_given,'MarkerFaceColor',[0.5,0.5,0.5],'Marker','o','Color',[0.5,0.5,0.5],'LineWidth',2,"LineStyle","none");

% % uninfected too, 5H community
% plot((1:5)+gap,gr_community,'o','MarkerSize',market_size_given,'MarkerFaceColor',color_green,'Color',color_green,'LineWidth',2);
% %errorbar((1:5),r_onestep,r_onestep_error)
% plot((1:5)-gap,gr_inf_community,'o','MarkerSize',market_size_given,'MarkerFaceColor','k','Color','k','LineWidth',2);
% %plot((1:5), [0.18, 0.22, 0.29, 0.67, 0.53], 'o','MarkerSize',market_size_given,'MarkerFaceColor','b','Color','b','LineWidth',2);

errorbar((1:5)+gap,gr_community,std_community,'o','MarkerSize',market_size_given,'MarkerFaceColor',color_green,'Color',color_green,'LineWidth',2);
errorbar((1:5)-gap,gr_inf_community,std_inf_community,'o','MarkerSize',market_size_given,'MarkerFaceColor','k','Color','k','LineWidth',2);


set(gca,'FontSize',20,'FontName','times');
xticks(1:5);
xticklabels({'CBA 4','CBA 18','CBA 38','PSA H100','PSA 13-15'});
set(gca,'XTickLabelRotation',90)
ylabel({'Growth rates','r (cells/hr)'});
xline(3.5,'--k',LineWidth=2.5);
title("Exponential fits to first 4 hours of the community")

legend("Phage-free control","Phage-infected community");
box on;