clc;
clear all;

red = [110, 27, 7]/255;
orange = [222, 133, 18]/255;
yellow = [194, 176, 19]/255;
green = [32, 117, 13]/255;
blue = [14, 97, 117]/255;

%% PSA H100;
load('./../../community/inhibition/H100.mat');

close all;
figure(1)
subplot(2,5,1)
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);

errorbar(time,mean(HP1_R),std(HP1_R),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(HP1_C),std(HP1_C),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('PSA H100 -- PSA HP1');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);




subplot(2,5,2)
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(HS6_R),std(HS6_R),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(HS6_C),std(HS6_C),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('PSA H100 -- PSA HS6');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);


subplot(2,5,3)
load('./../../community/inhibition/13-15.mat');
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(HP1_R),std(HP1_R),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(HP1_C),std(HP1_C),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('PSA 13-15 -- PSA HP1');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);


subplot(2,5,4)
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(HS6_R),std(HS6_R),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(HS6_C),std(HS6_C),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('PSA 13-15 -- PSA HS6');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);



subplot(2,5,6)
load('./../../community/inhibition/CBA4.mat')
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(R_18_3),std(R_18_3),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(C_18_3),std(C_18_3),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('CBA4 -- \phi18:3');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);


subplot(2,5,7)
load('./../../community/inhibition/CBA18.mat')
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(R_18_2),std(R_18_2),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(C_18_2),std(C_18_2),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('CBA18 -- \phi18:2');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);




subplot(2,5,8)
load('./../../community/inhibition/CBA18.mat')
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(R_18_3),std(R_18_3),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(C_18_3),std(C_18_3),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('CBA18 -- \phi18:3');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);




subplot(2,5,9)

errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(R_38_1),std(R_38_1),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(C_38_1),std(C_38_1),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('CBA18 -- \phi38:1');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);



subplot(2,5,10)
load('./../../community/inhibition/CBA38.mat')
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(R_38_1),std(R_38_1),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(C_38_1),std(C_38_1),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);

legend('control','infection from mock community','no phage control from mock community','added to pairwise infection','added to no phage control')

xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('CBA38 -- \phi38:1');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);






%% 


load('./../../community/inhibition/T12/H100_t12.mat');


figure(2)
subplot(2,5,1)
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(HP1_R),std(HP1_R),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(HP1_C),std(HP1_C),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('PSA H100 -- PSA HP1');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);




subplot(2,5,2)
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(HS6_R),std(HS6_R),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(HS6_C),std(HS6_C),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('PSA H100 -- PSA HS6');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);


subplot(2,5,3)
load('./../../community/inhibition/T12/13-15_t12.mat');
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(HP1_R),std(HP1_R),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(HP1_C),std(HP1_C),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('PSA 13-15 -- PSA HP1');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);


subplot(2,5,4)
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(HS6_R),std(HS6_R),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(HS6_C),std(HS6_C),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('PSA 13-15 -- PSA HS6');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);



subplot(2,5,6)
load('./../../community/inhibition/T12/CBA4_t12.mat')
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(R_18_3),std(R_18_3),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(C_18_3),std(C_18_3),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('CBA4 -- \phi18:3');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);


subplot(2,5,7)
load('./../../community/inhibition/T12/CBA18_t12.mat')
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(R_18_2),std(R_18_2),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(C_18_2),std(C_18_2),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('CBA18 -- \phi18:2');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);




subplot(2,5,8)
load('./../../community/inhibition/T12/CBA18_t12.mat')
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(R_18_3),std(R_18_3),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(C_18_3),std(C_18_3),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('CBA18 -- \phi18:3');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);




subplot(2,5,9)

errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(R_38_1),std(R_38_1),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(C_38_1),std(C_38_1),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('CBA18 -- \phi38:1');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);



subplot(2,5,10)
load('./../../community/inhibition/T12/CBA38_t12.mat')
errorbar(time,mean(control),std(control),'Marker','square', ...
    'MarkerEdgeColor',red,'MarkerFaceColor',red, ...
    'MarkerSize',10,LineWidth=2',Color=red);

hold on;

errorbar(time,mean(R),std(R),'Marker','square', ...
    'MarkerEdgeColor',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',10,LineWidth=2',Color=orange);

errorbar(time,mean(C),std(C),'Marker','square', ...
    'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow, ...
    'MarkerSize',10,LineWidth=2',Color=yellow);


errorbar(time,mean(R_38_1),std(R_38_1),'Marker','^', ...
    'MarkerEdgeColor',green,'MarkerFaceColor',green, ...
    'MarkerSize',10,LineWidth=2',Color=green);

errorbar(time,mean(C_38_1),std(C_38_1),'Marker','^', ...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',10,LineWidth=2',Color=blue);

legend('control','infection from mock community','no phage control from mock community','added to pairwise infection','added to no phage control')


xlabel('Time (hr)');
ylabel('OD_{600 nm}');
title('CBA38 -- \phi38:1');
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);










