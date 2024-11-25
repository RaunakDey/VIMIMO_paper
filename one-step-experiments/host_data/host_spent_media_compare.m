clc;
clear all;

% will think about the order later.

blank_correction = 1; % put 0 to turn off
blank = [0	-0.001	0.0002	0.0006	0.0006	0	0.0004	-0.0006	0	-0.0008	-0.0002	-0.0004	0.0026	-0.0002	0.003	0.0024	0.0016	0.0014	0.002	0.0024	0.0012	0.0028	0.0034	0.0022	0.0016	0.0046	0.0028	0.003];
time2 = [0	35	70	105	140	175	210	245	280	315	350	385	420	455	490	525	560	595	630	665	700	735	770	805	840	875	910	945]/60;
time = 0:19;

% we assume that both the experiments have same conversion units.
conversion_unit = [2.83e9/0.0577, 1.39e9/0.05055, 2.47e+09/0.06525, 1.77E+09/0.13815, 2.55E+09/0.10615];


%% cba 4

load("cba4.mat");
subplot(2,5,1);
plot(time2, conversion_unit(1)*( no_phage - blank_correction*blank) ,'ko','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
plot(time2, (phage_18_3 - blank_correction*blank)*conversion_unit(1),'mo', 'MarkerEdgeColor','m','MarkerFaceColor','m');



set(gca,'Yscale','log');
set(gca,'FontName','Times');
set(gca,'FontSize',20);
xlabel('Time (hrs)');
ylabel('Host density (OD)');
%ylim([1e-4, 5e-1]);
title('18:3 on CBA 4');
xticks([0 5 10 15]);
%yticks([1e-4 1e-3 1e-2 1e-1]);
ylim([1e6 1e10]);
yticks([1e6 1e7 1e8 1e9 1e10]);
ylabel('Host density (CFU/ml)')

load("./host_spent_media/CBA4-18-3.mat");
plot(time,conversion_unit(1)*debris_added,'go', 'MarkerEdgeColor','g','MarkerFaceColor','g');
xlim([0 15.75]);


%% CBA 18

load("cba18.mat");

subplot(2,5,2);
plot(time2,(no_phage - blank_correction*blank)*conversion_unit(2),'ko','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
plot(time2, (phage_18_2 - blank_correction*blank)*conversion_unit(2),'mo', 'MarkerEdgeColor','m','MarkerFaceColor','m');
set(gca,'Yscale','log');
set(gca,'FontName','Times');
set(gca,'FontSize',20);
xlabel('Time (hrs)');
ylabel('Host density (OD)');
%ylim([1e-4, 5e-1]);
title('18:2 on CBA 18');
xticks([0 5 10 15]);
%yticks([1e-4 1e-3 1e-2 1e-1]);
ylim([1e6 1e10]);
yticks([1e6 1e7 1e8 1e9 1e10]);
ylabel('Host density (CFU/ml)')
load("./host_spent_media/CBA18-18-2.mat");
plot(time,conversion_unit(2)*debris_added,'go', 'MarkerEdgeColor','g','MarkerFaceColor','g');
xlim([0 15.75]);


subplot(2,5,3);
plot(time2,(no_phage- blank_correction*blank)*conversion_unit(2),'ko','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
plot(time2, (phage_18_3- blank_correction*blank)*conversion_unit(2),'mo', 'MarkerEdgeColor','m','MarkerFaceColor','m');
set(gca,'Yscale','log');
set(gca,'FontName','Times');
set(gca,'FontSize',20);
xlabel('Time (hrs)');
ylabel('Host density (OD)');
%ylim([1e-4, 5e-1]);
title('18:3 on CBA 18');
xticks([0 5 10 15]);
%yticks([1e-4 1e-3 1e-2 1e-1]);
ylim([1e6 1e10]);
yticks([1e6 1e7 1e8 1e9 1e10]);
ylabel('Host density (CFU/ml)')
load("./host_spent_media/CBA18-18-3.mat");
plot(time,conversion_unit(2)*debris_added,'go', 'MarkerEdgeColor','g','MarkerFaceColor','g');
xlim([0 15.75]);


subplot(2,5,4);
plot(time2,(no_phage- blank_correction*blank)*conversion_unit(2),'ko','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
plot(time2, (phage_38_1- blank_correction*blank)*conversion_unit(2),'mo', 'MarkerEdgeColor','m','MarkerFaceColor','m');
set(gca,'Yscale','log');
set(gca,'FontName','Times');
set(gca,'FontSize',20);
xlabel('Time (hrs)');
ylabel('Host density (OD)');
%ylim([1e-4, 5e-1]);
title('38:1 on CBA 18');
xticks([0 5 10 15]);
%yticks([1e-4 1e-3 1e-2 1e-1]);
ylim([1e6 1e10]);
yticks([1e6 1e7 1e8 1e9 1e10]);
ylabel('Host density (CFU/ml)')
load("./host_spent_media/CBA18-38-1.mat")
plot(time,conversion_unit(2)*debris_added,'go', 'MarkerEdgeColor','g','MarkerFaceColor','g');
xlim([0 15.75]);

%% cba 38

load("cba38.mat");

subplot(2,5,5)
plot(time2,(no_phage- blank_correction*blank)*conversion_unit(3),'ko','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
plot(time2,(phage_38_1 - blank_correction*blank)*conversion_unit(3),'mo', 'MarkerEdgeColor','m','MarkerFaceColor','m');
set(gca,'Yscale','log');
set(gca,'FontName','Times');
set(gca,'FontSize',20);
xlabel('Time (hrs)');
ylabel('Host density (OD)');
%ylim([1e-4, 5e-1]);
title('38:1 on CBA 38');
xticks([0 5 10 15]);
%yticks([1e-4 1e-3 1e-2 1e-1]);
ylim([1e6 1e10]);
yticks([1e6 1e7 1e8 1e9 1e10]);
ylabel('Host density (CFU/ml)')
load("./host_spent_media/CBA38-38-1.mat")
plot(time,conversion_unit(3)*debris_added,'go', 'MarkerEdgeColor','g','MarkerFaceColor','g');
xlim([0 15.75]);


%% psa h100

load("psa_h100.mat");

subplot(2,5,6)
plot(time2,(no_phage - blank_correction*blank)*conversion_unit(4),'ko','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
plot(time2,(phage_hp1 - blank_correction*blank)*conversion_unit(4),'mo', 'MarkerEdgeColor','m','MarkerFaceColor','m');
set(gca,'Yscale','log');
set(gca,'FontName','Times');
set(gca,'FontSize',20);
xlabel('Time (hrs)');
ylabel('Host density (OD)');
%ylim([1e-4, 5e-1]);
title('HP1 on PSA H100');
xticks([0 5 10 15]);
%yticks([1e-4 1e-3 1e-2 1e-1]);
ylim([1e6 1e10]);
yticks([1e6 1e7 1e8 1e9 1e10]);
ylabel('Host density (CFU/ml)')
load("./host_spent_media/h100_hp1.mat")
plot(time,conversion_unit(4)*debris_added,'go', 'MarkerEdgeColor','g','MarkerFaceColor','g');
xlim([0 15.75]);


subplot(2,5,7)
plot(time2,(no_phage - blank_correction*blank)*conversion_unit(4),'ko','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
plot(time2,(phage_hs6 - blank_correction*blank)*conversion_unit(4),'mo', 'MarkerEdgeColor','m','MarkerFaceColor','m');
set(gca,'Yscale','log');
set(gca,'FontName','Times');
set(gca,'FontSize',20);
xlabel('Time (hrs)');
ylabel('Host density (OD)');
%ylim([1e-4, 5e-1]);
title('HS6 on PSA H100');
xticks([0 5 10 15]);
%yticks([1e-4 1e-3 1e-2 1e-1]);
ylim([1e6 1e10]);
yticks([1e6 1e7 1e8 1e9 1e10]);
ylabel('Host density (CFU/ml)')
load("./host_spent_media/h100_hs6.mat")
plot(time,conversion_unit(4)*debris_added,'go', 'MarkerEdgeColor','g','MarkerFaceColor','g');
xlim([0 15.75]);

%% psa 13-15


load("psa_13_15.mat");

subplot(2,5,8)
plot(time2,(no_phage - blank_correction*blank)*conversion_unit(5),'ko','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
plot(time2,(phage_hp1 - blank_correction*blank)*conversion_unit(5),'mo', 'MarkerEdgeColor','m','MarkerFaceColor','m');
set(gca,'Yscale','log');
set(gca,'FontName','Times');
set(gca,'FontSize',20);
xlabel('Time (hrs)');
ylabel('Host density (OD)');
%ylim([1e-4, 5e-1]);
title('HP1 on PSA 13-15');
xticks([0 5 10 15]);
%yticks([1e-4 1e-3 1e-2 1e-1]);
ylim([1e6 1e10]);
yticks([1e6 1e7 1e8 1e9 1e10]);
ylabel('Host density (CFU/ml)')
load("./host_spent_media/1315_hp1.mat")
plot(time,conversion_unit(5)*debris_added,'go', 'MarkerEdgeColor','g','MarkerFaceColor','g');
xlim([0 15.75]);


subplot(2,5,9)




plot(time2,(no_phage(3,:) - blank_correction*blank)*conversion_unit(5),'ko','MarkerEdgeColor','k','MarkerFaceColor','k','DisplayName', 'phage free single host');
hold on;
plot(time2,(phage_hs6(3,:) - blank_correction*blank)*conversion_unit(5),'mo','MarkerEdgeColor','m','MarkerFaceColor','m','DisplayName', 'single phage infected single hosts');

load("./host_spent_media/1315_hs6.mat")
plot(time,conversion_unit(2)*debris_added(3,:),'go', 'MarkerEdgeColor','g','MarkerFaceColor','g','DisplayName', 'single phage infected single hosts + spent media from phage infected 5h5v community at 12 hrs');
xlim([0 15.75]);


plot(time2,(no_phage(1:2,:) - blank_correction*blank)*conversion_unit(5),'ko','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
plot(time2,(phage_hs6(1:2,:) - blank_correction*blank)*conversion_unit(5),'mo','MarkerEdgeColor','m','MarkerFaceColor','m');
plot(time,conversion_unit(5)*debris_added(1:2,:),'go', 'MarkerEdgeColor','g','MarkerFaceColor','g')

set(gca,'Yscale','log');
set(gca,'FontName','Times');
set(gca,'FontSize',20);
xlabel('Time (hrs)');
ylabel('Host density (OD)');

ylabel('Host density (CFU/ml)')
%ylim(conversion_unit(5)*[1e-4, 5e-1]);
title('HS6 on PSA 13-15');
xticks([0 5 10 15]);
%yticks([1e-4 1e-3 1e-2 1e-1]);
ylim([1e6 1e10]);
yticks([1e6 1e7 1e8 1e9 1e10]);



% Display the legend for this subplot
legend('pphage free single host','single phage infected single hosts','single phage infected single hosts spent media from 12 hrs community');
legend('show');