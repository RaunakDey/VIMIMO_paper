clear all;
clc;


% Please not fig 2a is a schematics not a plot.

load('./../community/data/triplicate_data.mat');
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

%% plots
figure(1)
tiledlayout(2,1,"TileSpacing","none");

% Tile 1
nexttile
errorbar(time, mean(host1'), std(host1'),Marker="o",MarkerFaceColor =blue1,MarkerSize=10,Color=blue1,LineStyle="-",LineWidth=2); hold on;
errorbar(time, mean(host2'), std(host2'),Marker="o",MarkerFaceColor =blue2,MarkerSize=10,Color=blue2,LineStyle="-",LineWidth=2); 
errorbar(time, mean(host3'), std(host3'),Marker="o",MarkerFaceColor =blue3,MarkerSize=10,Color=blue3,LineStyle="-",LineWidth=2); 
errorbar(time, mean(host4'), std(host4'),Marker="o",MarkerFaceColor =blue4,MarkerSize=10,Color=blue4,LineStyle="-",LineWidth=2); 
errorbar(time, mean(host5'), std(host5'),Marker="o",MarkerFaceColor =blue5,MarkerSize=10,Color=blue5,LineStyle="-",LineWidth=2); 

hold on;
set(gca,'FontName','Times');
set(gca,'FontSize',30);
ylabel('Bacterial density (cells/ml)');
set(gca,'YScale','log');
ylim([1e5 1e8]);
yticks([ 1e6, 1e7, 1e8]);
xticks([]);
legend('CBA4','CBA18', 'CBA38','PSA H100','PSA 13-15')


% Tile 2
nexttile
errorbar(time, mean(virus1'), std(virus1'),Marker="o",MarkerFaceColor =red1, MarkerSize=10, Color=red1,LineStyle="-",LineWidth=2); hold on;
errorbar(time, mean(virus2'), std(virus2'),Marker="o",MarkerFaceColor =red2, MarkerSize=10, Color=red2,LineStyle="-",LineWidth=2); 
errorbar(time, mean(virus3'), std(virus3'),Marker="o",MarkerFaceColor =red3, MarkerSize=10, Color=red3,LineStyle="-",LineWidth=2); 
errorbar(time, mean(virus4'), std(virus4'),Marker="o",MarkerFaceColor =red4, MarkerSize=10, Color=red4,LineStyle="-",LineWidth=2); 
errorbar(time, mean(virus5'), std(virus5'),Marker="o",MarkerFaceColor =red5, MarkerSize=10, Color=red5,LineStyle="-",LineWidth=2); 

hold on;
set(gca,'FontSize',30);
set(gca,'FontName','Times');
ylabel('Phage density (cells/ml)')
set(gca,'YScale','log');
xticks(0:1:16);
ylim([1e4 1e12]);
yticks([1e4, 1e5,1e6,1e7, 1e8, 1e9, 1e10, 1e11]);
legend('\phi18:2','\phi18:3', '\phi38:1','PSA HP1','PSA HS6')
xlabel('Time (hrs)')

%%

default_msize = 10;

figure(2)
subplot(2,5,1)
plot(time,host1(:,1),Color=blue1, MarkerEdgeColor='k',MarkerFaceColor=blue1, LineStyle="none",Marker="o",MarkerSize=default_msize );
hold on;
plot(time,host1(:,2),Color=blue1 ,MarkerEdgeColor='k',MarkerFaceColor=blue1,LineStyle="none",Marker="o",MarkerSize=default_msize );
plot(time,host1(:,3),Color=blue1,MarkerEdgeColor='k',MarkerFaceColor=blue1, LineStyle="none",Marker="o",MarkerSize=default_msize );
set(gca,'FontSize',20);
set(gca,'FontName','Times');
xlabel('Time (hrs)');
ylabel('Bacterial density (cells/ml)');
ylim([1e5 1e8]);
yticks([1e5, 1e6,1e7,1e8]);
set(gca, 'YScale','log');
xlim([0 16]);
xticks(0:2:16);
title('CBA 4');
xtickangle(90);



subplot(2,5,2)
plot(time,host2(:,1),Color=blue2, MarkerEdgeColor='k',MarkerFaceColor=blue2, LineStyle="none",Marker="o",MarkerSize=default_msize );
hold on;
plot(time,host2(:,2),Color=blue2 ,MarkerEdgeColor='k',MarkerFaceColor=blue2,LineStyle="none",Marker="o",MarkerSize=default_msize );
plot(time,host2(:,3),Color=blue2,MarkerEdgeColor='k',MarkerFaceColor=blue2, LineStyle="none",Marker="o",MarkerSize=default_msize );
set(gca,'FontSize',20);
set(gca,'FontName','Times');
xlabel('Time (hrs)');
ylim([1e5 1e8]);
yticks([1e5, 1e6,1e7,1e8]);
set(gca, 'YScale','log');
xlim([0 16]);
xticks(0:2:16);
title('CBA 18');
xtickangle(90);



subplot(2,5,3)
plot(time,host3(:,1),Color=blue3, MarkerEdgeColor='k',MarkerFaceColor=blue3, LineStyle="none",Marker="o",MarkerSize=default_msize );
hold on;
plot(time,host3(:,2),Color=blue3 ,MarkerEdgeColor='k',MarkerFaceColor=blue3,LineStyle="none",Marker="o",MarkerSize=default_msize );
plot(time,host3(:,3),Color=blue3,MarkerEdgeColor='k',MarkerFaceColor=blue3, LineStyle="none",Marker="o",MarkerSize=default_msize );
set(gca,'FontSize',20);
set(gca,'FontName','Times');
xlabel('Time (hrs)');
ylim([1e5 1e8]);
yticks([1e5, 1e6,1e7,1e8]);
set(gca, 'YScale','log');
xlim([0 16]);
xticks(0:2:16);
title('CBA 38');
xtickangle(90);



subplot(2,5,4)
plot(time,host4(:,1),Color=blue4, MarkerEdgeColor='k',MarkerFaceColor=blue4, LineStyle="none",Marker="o",MarkerSize=default_msize );
hold on;
plot(time,host4(:,2),Color=blue4 ,MarkerEdgeColor='k',MarkerFaceColor=blue4,LineStyle="none",Marker="o",MarkerSize=default_msize );
plot(time,host4(:,3),Color=blue4,MarkerEdgeColor='k',MarkerFaceColor=blue4, LineStyle="none",Marker="o",MarkerSize=default_msize );
set(gca,'FontSize',20);
set(gca,'FontName','Times');
xlabel('Time (hrs)');
ylim([1e5 1e8]);
yticks([1e5, 1e6,1e7,1e8]);
set(gca, 'YScale','log');
xlim([0 16]);
xticks(0:2:16);
title('PSA H100');
xtickangle(90);



subplot(2,5,5)
plot(time,host5(:,1),Color=blue5, MarkerEdgeColor='k',MarkerFaceColor=blue5, LineStyle="none",Marker="o",MarkerSize=default_msize );
hold on;
plot(time,host5(:,2),Color=blue5 ,MarkerEdgeColor='k',MarkerFaceColor=blue5,LineStyle="none",Marker="o",MarkerSize=default_msize );
plot(time,host5(:,3),Color=blue5,MarkerEdgeColor='k',MarkerFaceColor=blue5, LineStyle="none",Marker="o",MarkerSize=default_msize );
set(gca,'FontSize',20);
set(gca,'FontName','Times');
xlabel('Time (hrs)');
ylim([1e5 1e8]);
yticks([1e5, 1e6,1e7,1e8]);
set(gca, 'YScale','log');
xlim([0 16]);
xticks(0:2:16);
title('PSA 13-15');
xtickangle(90);




subplot(2,5,6)
plot(time,virus1(:,1),Color=red1, MarkerEdgeColor='k',MarkerFaceColor=red1, LineStyle="none",Marker="o",MarkerSize=default_msize );
hold on;
plot(time,virus1(:,2),Color=red1 ,MarkerEdgeColor='k',MarkerFaceColor=red1,LineStyle="none",Marker="o",MarkerSize=default_msize );
plot(time,virus1(:,3),Color=red1,MarkerEdgeColor='k',MarkerFaceColor=red1, LineStyle="none",Marker="o",MarkerSize=default_msize );
set(gca,'FontSize',20);
set(gca,'FontName','Times');
xlabel('Time (hrs)');
ylabel('Phage density (virions/ml)');
ylim([1e4 1e11]);
yticks([1e4, 1e6,1e8,1e10]);
set(gca, 'YScale','log');
xlim([0 16]);
xticks(0:2:16);
title('\phi18:2');
xtickangle(90);



subplot(2,5,7)
plot(time,virus2(:,1),Color=red2, MarkerEdgeColor='k',MarkerFaceColor=red2, LineStyle="none",Marker="o",MarkerSize=default_msize );
hold on;
plot(time,virus2(:,2),Color=red2 ,MarkerEdgeColor='k',MarkerFaceColor=red2,LineStyle="none",Marker="o",MarkerSize=default_msize );
plot(time,virus2(:,3),Color=red2,MarkerEdgeColor='k',MarkerFaceColor=red2, LineStyle="none",Marker="o",MarkerSize=default_msize );
set(gca,'FontSize',20);
set(gca,'FontName','Times');
xlabel('Time (hrs)');
ylim([1e4 1e11]);
yticks([1e4, 1e6,1e8,1e10]);
set(gca, 'YScale','log');
xlim([0 16]);
xticks(0:2:16);
title('\phi18:3');
xtickangle(90);


subplot(2,5,8)
plot(time,virus3(:,1),Color=red3, MarkerEdgeColor='k',MarkerFaceColor=red3, LineStyle="none",Marker="o",MarkerSize=default_msize );
hold on;
plot(time,virus3(:,2),Color=red3 ,MarkerEdgeColor='k',MarkerFaceColor=red3,LineStyle="none",Marker="o",MarkerSize=default_msize );
plot(time,virus3(:,3),Color=red3,MarkerEdgeColor='k',MarkerFaceColor=red3, LineStyle="none",Marker="o",MarkerSize=default_msize );
set(gca,'FontSize',20);
set(gca,'FontName','Times');
xlabel('Time (hrs)');
ylim([1e4 1e11]);
yticks([1e4, 1e6,1e8,1e10]);
set(gca, 'YScale','log');
xlim([0 16]);
xticks(0:2:16);
title('\phi38:1');
xtickangle(90);


subplot(2,5,9)
plot(time,virus4(:,1),Color=red4, MarkerEdgeColor='k',MarkerFaceColor=red4, LineStyle="none",Marker="o",MarkerSize=default_msize );
hold on;
plot(time,virus4(:,2),Color=red4 ,MarkerEdgeColor='k',MarkerFaceColor=red4,LineStyle="none",Marker="o",MarkerSize=default_msize );
plot(time,virus4(:,3),Color=red4,MarkerEdgeColor='k',MarkerFaceColor=red4, LineStyle="none",Marker="o",MarkerSize=default_msize );
set(gca,'FontSize',20);
set(gca,'FontName','Times');
xlabel('Time (hrs)');
ylim([1e4 1e11]);
yticks([1e4, 1e6,1e8,1e10]);
set(gca, 'YScale','log');
xlim([0 16]);
xticks(0:2:16);
title('PSA HP1');
xtickangle(90);


subplot(2,5,10)
plot(time,virus5(:,1),Color=red5, MarkerEdgeColor='k',MarkerFaceColor=red5, LineStyle="none",Marker="o",MarkerSize=default_msize );
hold on;
plot(time,virus5(:,2),Color=red5 ,MarkerEdgeColor='k',MarkerFaceColor=red5,LineStyle="none",Marker="o",MarkerSize=default_msize );
plot(time,virus5(:,3),Color=red5,MarkerEdgeColor='k',MarkerFaceColor=red5, LineStyle="none",Marker="o",MarkerSize=default_msize );
set(gca,'FontSize',20);
set(gca,'FontName','Times');
xlabel('Time (hrs)');
ylim([1e4 1e11]);
yticks([1e4, 1e6,1e8,1e10]);
set(gca, 'YScale','log');
xlim([0 16]);
xticks(0:2:16);
title('PSA HS6');
xtickangle(90);


