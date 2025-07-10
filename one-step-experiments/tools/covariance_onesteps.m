clc;
clear all;
close all;

%% load data

%load('./results_same_phi/CBA18-3_4_10-inferred.mat');
load('./results_same_phi/HP1_13-15_3-inferred.mat');


transient_id = 3000;

figure;
subplot(1,3,1)
i= 2;
j= 3;
pearson = corrcoef(chain(transient_id:end,i),chain(transient_id:end,j));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

plot(chain(transient_id:end,j),chain(transient_id:end,i),'.','Color',color_input);
set(gca,'FontSize',14);
axis square;
ylabel('\phi (mL/hr)');
xlabel('\tau (hr)');



subplot(1,3,2)
i= 3;
j= 4;
pearson = corrcoef(chain(transient_id:end,i),chain(transient_id:end,j));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

plot(chain(transient_id:end,j),chain(transient_id:end,i),'.','Color',color_input);
set(gca,'FontSize',14);
axis square;

ylabel('\tau (hr)');
xlabel('\beta (virions/cells)');


subplot(1,3,3)
i= 4;
j= 2;
pearson = corrcoef(chain(transient_id:end,i),chain(transient_id:end,j));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

plot(chain(transient_id:end,j),chain(transient_id:end,i),'.','Color',color_input);
set(gca,'FontSize',14);
axis square;
ylabel('\tau (hr)');
xlabel('\phi (mL/hr)');
