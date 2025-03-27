clear all;
clc;

time = 0:19;
green = [32, 117, 13]/255;

%% load

figure(1)


load("h100_hp1.mat");
subplot(2,5,1)
errorbar(time, mean(only_virus),std(only_virus),'k','LineWidth',2,Marker='o',MarkerEdgeColor='k',MarkerFaceColor='k',MarkerSize=8);
hold on;
errorbar(time, mean(debris_added),std(debris_added),'Color',green,'LineWidth',2,Marker='square',MarkerEdgeColor=green,MarkerFaceColor=green,MarkerSize=8);
set(gca,'YScale','log');
set(gca, 'FontName','times')
set(gca, 'FontSize',20)
title('H100 -- HP1')
clear only_virus debris_added;
ylim([1e-4 5e-1]);
yticks([1e-4 1e-3 1e-2 1e-1]);
xlabel('Time (hrs)');
ylabel('Host density (OD_{600nm})');



load("h100_hs6.mat");
subplot(2,5,2)
errorbar(time, mean(only_virus),std(only_virus),'k','LineWidth',2,Marker='o',MarkerEdgeColor='k',MarkerFaceColor='k',MarkerSize=8);
hold on;
errorbar(time, mean(debris_added),std(debris_added),'Color',green,'LineWidth',2,Marker='square',MarkerEdgeColor=green,MarkerFaceColor=green,MarkerSize=8);
set(gca,'YScale','log');
set(gca, 'FontName','times')
set(gca, 'FontSize',20)
title('H100 -- HS6')
clear only_virus debris_added;
ylim([1e-4 5e-1]);
yticks([1e-4 1e-3 1e-2 1e-1]);
xlabel('Time (hrs)');
ylabel('Host density (OD_{600nm})');




load("1315_hp1.mat"); %change here
subplot(2,5,3) %change here
errorbar(time, mean(only_virus),std(only_virus),'k','LineWidth',2,Marker='o',MarkerEdgeColor='k',MarkerFaceColor='k',MarkerSize=8);
hold on;
errorbar(time, mean(debris_added),std(debris_added),'Color',green,'LineWidth',2,Marker='square',MarkerEdgeColor=green,MarkerFaceColor=green,MarkerSize=8);
set(gca,'YScale','log');
set(gca, 'FontName','times')
set(gca, 'FontSize',20)
title('H13-15 -- HP1') %change here.
clear only_virus debris_added;
ylim([1e-4 5e-1]);
yticks([1e-4 1e-3 1e-2 1e-1]);
xlabel('Time (hrs)');
ylabel('Host density (OD_{600nm})');



load("1315_hs6.mat"); %change here
subplot(2,5,4) %change here
errorbar(time, mean(only_virus),std(only_virus),'k','LineWidth',2,Marker='o',MarkerEdgeColor='k',MarkerFaceColor='k',MarkerSize=8);
hold on;
errorbar(time, mean(debris_added),std(debris_added),'Color',green,'LineWidth',2,Marker='square',MarkerEdgeColor=green,MarkerFaceColor=green,MarkerSize=8);
set(gca,'YScale','log');
set(gca, 'FontName','times')
set(gca, 'FontSize',20)
title('H13-15 -- HS6') %change here.
clear only_virus debris_added;
ylim([1e-4 5e-1]);
yticks([1e-4 1e-3 1e-2 1e-1]);
xlabel('Time (hrs)');
ylabel('Host density (OD_{600nm})');



load("CBA4-18-3.mat"); %change here
subplot(2,5,6) %change here
errorbar(time, mean(only_virus),std(only_virus),'k','LineWidth',2,Marker='o',MarkerEdgeColor='k',MarkerFaceColor='k',MarkerSize=8);
hold on;
errorbar(time, mean(debris_added),std(debris_added),'Color',green,'LineWidth',2,Marker='square',MarkerEdgeColor=green,MarkerFaceColor=green,MarkerSize=8);
set(gca,'YScale','log');
set(gca, 'FontName','times')
set(gca, 'FontSize',20)
title('CBA4 -- CBA 18:3') %change here.
clear only_virus debris_added;
ylim([1e-4 5e-1]);
yticks([1e-4 1e-3 1e-2 1e-1]);
xlabel('Time (hrs)');
ylabel('Host density (OD_{600nm})');



load("CBA18-18-2.mat"); %change here
subplot(2,5,7) %change here
errorbar(time, mean(only_virus),std(only_virus),'k','LineWidth',2,Marker='o',MarkerEdgeColor='k',MarkerFaceColor='k',MarkerSize=8);
hold on;
errorbar(time, mean(debris_added),std(debris_added),'Color',green,'LineWidth',2,Marker='square',MarkerEdgeColor=green,MarkerFaceColor=green,MarkerSize=8);
set(gca,'YScale','log');
set(gca, 'FontName','times')
set(gca, 'FontSize',20)
title('CBA18 -- CBA 18:2') %change here.
clear only_virus debris_added;
ylim([1e-4 5e-1]);
yticks([1e-4 1e-3 1e-2 1e-1]);
xlabel('Time (hrs)');
ylabel('Host density (OD_{600nm})');




load("CBA18-18-3.mat"); %change here
subplot(2,5,8) %change here
errorbar(time, mean(only_virus),std(only_virus),'k','LineWidth',2,Marker='o',MarkerEdgeColor='k',MarkerFaceColor='k',MarkerSize=8);
hold on;
errorbar(time, mean(debris_added),std(debris_added),'Color',green,'LineWidth',2,Marker='square',MarkerEdgeColor=green,MarkerFaceColor=green,MarkerSize=8);
set(gca,'YScale','log');
set(gca, 'FontName','times')
set(gca, 'FontSize',20)
title('CBA18 -- CBA 18:3') %change here.
clear only_virus debris_added;
ylim([1e-4 5e-1]);
yticks([1e-4 1e-3 1e-2 1e-1]);
xlabel('Time (hrs)');
ylabel('Host density (OD_{600nm})');




load("CBA18-38-1.mat"); %change here
subplot(2,5,9) %change here
errorbar(time, mean(only_virus),std(only_virus),'k','LineWidth',2,Marker='o',MarkerEdgeColor='k',MarkerFaceColor='k',MarkerSize=8);
hold on;
errorbar(time, mean(debris_added),std(debris_added),'Color',green,'LineWidth',2,Marker='square',MarkerEdgeColor=green,MarkerFaceColor=green,MarkerSize=8);
set(gca,'YScale','log');
set(gca, 'FontName','times')
set(gca, 'FontSize',20)
title('CBA18 -- CBA 38:1') %change here.
clear only_virus debris_added;
ylim([1e-4 5e-1]);
yticks([1e-4 1e-3 1e-2 1e-1]);
xlabel('Time (hrs)');
ylabel('Host density (OD_{600nm})');




load("CBA38-38-1.mat"); %change here
subplot(2,5,10) %change here
errorbar(time, mean(only_virus),std(only_virus),'k','LineWidth',2,Marker='o',MarkerEdgeColor='k',MarkerFaceColor='k',MarkerSize=8);
hold on;
errorbar(time, mean(debris_added),std(debris_added),'Color',green,'LineWidth',2,Marker='square',MarkerEdgeColor=green,MarkerFaceColor=green,MarkerSize=8);
set(gca,'YScale','log');
set(gca, 'FontName','times')
set(gca, 'FontSize',20)
title('CBA38 -- CBA 38:1') %change here.
clear only_virus debris_added;
ylim([1e-4 5e-1]);
yticks([1e-4 1e-3 1e-2 1e-1]);
xlabel('Time (hrs)');
ylabel('Host density (OD_{600nm})');


legend('Host + virus + spent media (T=0hr)','Host + virus + spent media (filtered/UV treated) (T=12hr)')