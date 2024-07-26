close all;
clear all;
clc;

%% load data

load('./../data/qpcr.mat');


load('./host_od.mat');
names = ['CBA4', 'CBA18', 'CBA38', 'PSA100', 'PSA13-15'];



od_total = sum(od_avg);

%%% calibration factor given by OSU experiments -- 2.97e9 OD 

calibration_factor = 2.97e9;

recalibrate = 0.3;
%recalibrate = 1;

cfu_total = calibration_factor * od_total;
cfu_total = cfu_total * recalibrate;

host_qpcr = data.ydata(:,1:5);
host_sum_with_phage = sum(host_qpcr');


%% figures

figure(1)
plot(time,cfu_total,'bo','MarkerSize',8,'MarkerFaceColor','b');
hold on;
plot(data.xdata, host_sum_with_phage,'r^','MarkerSize',8,'MarkerFaceColor','r')
xlabel('time (hr)');
ylabel('host density (cell/ml) ');
set(gca,'YScale','log');
set(gca,'FontSize',20);
legend('sum of host density without phage (CFU)','sum of host density with phages qpcr')
