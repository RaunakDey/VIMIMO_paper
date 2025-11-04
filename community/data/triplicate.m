clear all;
clc;

%% datasets.
load('./triplicate_data.mat');

%% plotting

host1 =host1*1e3;
host2 =host2*1e3;
host3 =host3*1e3;
host4 =host4*1e3;
host5 =host5*1e3;
virus1 = virus1*1e3;
virus2 = virus2*1e3;
virus3 = virus3*1e3;
virus4 = virus4*1e3;
virus5 = virus5*1e3;

time = time/60;




red1 = [1 0.1 0.1];
red2 = [0.8 0.3 0.3];
red3 = [0.6 0.4 0.4];

blue1 = [0.1 0.1 1];
blue2 = [0.4 0.4 0.8];

%% plots
figure(1)
tiledlayout(2,1,"TileSpacing","none");

% Tile 1
nexttile
errorbar(time, mean(host1'), std(host1'),Marker="o",MarkerFaceColor =red1,MarkerSize=10,Color=red1,LineStyle="-",LineWidth=2); hold on;
errorbar(time, mean(host2'), std(host2'),Marker="o",MarkerFaceColor =red2,MarkerSize=10,Color=red2,LineStyle="-",LineWidth=2); 
errorbar(time, mean(host3'), std(host3'),Marker="o",MarkerFaceColor =red3,MarkerSize=10,Color=red3,LineStyle="-",LineWidth=2); 
errorbar(time, mean(host4'), std(host4'),Marker="o",MarkerFaceColor =blue1,MarkerSize=10,Color=blue1,LineStyle="-",LineWidth=2); 
errorbar(time, mean(host5'), std(host5'),Marker="o",MarkerFaceColor =blue2,MarkerSize=10,Color=blue2,LineStyle="-",LineWidth=2); 

hold on;
set(gca,'FontName','Times');
set(gca,'FontSize',30);
ylabel('Host density (cells/ml)');
set(gca,'YScale','log');
ylim([1e5 1e8]);
yticks([ 1e6, 1e7, 1e8]);
xticks([]);
legend('CBA4','CBA18', 'CBA38','PSA H100','PSA 13-15')


% Tile 2
nexttile
errorbar(time, mean(virus1'), std(virus1'),Marker="*",MarkerFaceColor =red1, MarkerSize=10, Color=red1,LineStyle="-",LineWidth=2); hold on;
errorbar(time, mean(virus2'), std(virus2'),Marker="*",MarkerFaceColor =red2, MarkerSize=10,Color=red2,LineStyle="-",LineWidth=2); 
errorbar(time, mean(virus3'), std(virus3'),Marker="*",MarkerFaceColor =red3,MarkerSize=10,Color=red3,LineStyle="-",LineWidth=2); 
errorbar(time, mean(virus4'), std(virus4'),Marker="*",MarkerFaceColor =blue1,MarkerSize=10,Color=blue1,LineStyle="-",LineWidth=2); 
errorbar(time, mean(virus5'), std(virus5'),Marker="*",MarkerFaceColor =blue2,MarkerSize=10,Color=blue2,LineStyle="-",LineWidth=2); 

hold on;
set(gca,'FontSize',30);
set(gca,'FontName','Times');
ylabel('Phage density (cells/ml)')
set(gca,'YScale','log');
xticks(0:2:16);
ylim([1e4 1e12]);
yticks([1e4, 1e5,1e6,1e7, 1e8, 1e9, 1e10, 1e11]);
legend('\phi18:2','\phi18:3', '\phi38:1','PSA HP1','PSA HS6')
xlabel('Time (hrs)')


