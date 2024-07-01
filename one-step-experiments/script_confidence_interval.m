clear all;
clc;


%load("./results_same_phi/CBA18-2_18_5-inferred.mat");
%load("./results_same_phi/CBA18-3_4_10-inferred.mat");
load("./results_same_phi/CBA18-3_18_7-inferred.mat");
%load("./results_same_phi/CBA38-1_38_6-inferred.mat");
%load("./results_same_phi/HP1_H100_6-inferred.mat");
%load("./results_same_phi/HP1_13-15_3-inferred.mat");
%load("./results_same_phi/HS6_H100_3-inferred.mat");
%load("./results_same_phi/HS6_13-15_3-inferred.mat");

chain_effective = chain(5001:10000,:);

for i =1:3000
i
%NE = round(chain_effective(i,5));
y0(1) = S0;
y0(2:NE+2) = 0;
y0(NE+3) = V0;

[time_used,everything] = one_step_simulate_particular_points(time_free_phages,y0,chain_effective(i,1:5),NE,dilution_factor);
V(i,:) = everything(end,:);
clear y0;
end

Q = quantile(V,[0.4 0.75]);
max = Q(1,:);
min = Q(2,:);

save('./results_same_phi/'+name+'_'+seed+'-inferred.mat');


%% plotting

figure(5)
[~,something] = size(free_phages);


if something == 1
plot(time_free_phages,free_phages,'o','LineWidth', 2);
else
errorbar(time_free_phages,mean(free_phages,2),std(free_phages'),'LineWidth',2,'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',[70/255,130/255,180/255],Color=[70/255,130/255,180/255]);hold on

end
xlabel('time (hours)','interpreter','none')
ylabel("Free phages (" +string(cell2mat(labels.units(1,3)))+")",'interpreter','none');
title(string(labels.phage)+'-'+string(labels.host));
hold on;
plot(time,y_series_inferred(end,:),'LineWidth',2,Color=[171,193,157]./255);hold on;
    

time_2 = [time_used, fliplr(time_used)];
inBetween = [min, fliplr(max)];
fill(time_2, inBetween, [171,193,157]./255,'FaceAlpha',0.2,'LineStyle','none'); hold on;



legend('One step data','SEIV with inferred parameters','95% confidence interval','Location','northoutside');
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
str2 = {['Inferred'],['r = ',num2str(theta_inferred(1)),' cells/hr'], ['\phi = ',num2str(theta_inferred(2))], ['\tau = ',num2str(theta_inferred(3)),' hr'], ['\beta = ',num2str(theta_inferred(4))] ['N_E = ',num2str(round(theta_inferred(5)))]};
annotation('textbox', [0.7, 0.25, 0.1, 0.1], 'String', str2,'FontSize',11,'FitBoxToText','on');

% 
% 
% 
