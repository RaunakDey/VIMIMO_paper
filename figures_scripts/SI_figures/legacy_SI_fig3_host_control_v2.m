clear all;
clc;

%% datasets.

load('./../../community/data/nophage_control.mat');
time_no_phage = time;
load('./../../community/data/triplicate_data.mat');

%% plotting

host1 =host1*1e3;
host2 =host2*1e3;
host3 =host3*1e3;
host4 =host4*1e3;
host5 =host5*1e3;
CBA4 = CBA4*1e3;
CBA18 = CBA18*1e3;
CBA38 = CBA38*1e3;
H100 = H100*1e3;
PSA13_15 = PSA13_15*1e3;

time = time/60;


color_green = [171,193,157]./255;



blue1 = [179,205,224]/255;
blue2 = [100,151,177]/255;
blue3 = [0,91,150]/255;
blue4 = [3,10,208]/255;
blue5 = [0,0,75]/255;


%% fit growth rates 15.75 hours
g_low = zeros(1,5);
g_high = zeros(1,5);
growthrate = zeros(1,5);


bacteria = mean(CBA4,2);
% Fixed Initial Population (N0)
N0 = bacteria(1); % Set N0 to the first data point
% Logistic Growth Function with Fixed N0
logisticFunFixedN0 = @(params, t) params(1) ./ (1 + ((params(1) - N0) / N0) .* exp(-params(2) * t));
% Initial Guess for Parameters: [Carrying capacity K, Growth rate r]
initialParams = [max(bacteria), 0.1];
% Fit the model using non-linear regression with lsqcurvefit
[fitParams, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(logisticFunFixedN0, initialParams, time_no_phage, bacteria);
% Extract fitted parameters
K = fitParams(1);      % Carrying capacity
r = fitParams(2);      % Growth rate
% Calculate standard errors and confidence intervals
ci = nlparci(fitParams, residual, 'jacobian', jacobian); % 95% confidence intervals

growthrate(1) = r;
g_low(1) = ci(2,1);
g_high(1) = ci(2,2);
fittedBacteria(:,1) = logisticFunFixedN0(fitParams, time);

% Display Results
fprintf('Fitted Parameters with Confidence Intervals:\n');
fprintf('Carrying Capacity (K): %.3f (95%% CI: [%.3f, %.3f])\n', K, ci(1,1), ci(1,2));
fprintf('Growth Rate (r): %.3f (95%% CI: [%.3f, %.3f])\n', r, ci(2,1), ci(2,2));

% Plot the results
figure;
hold on;
scatter(time, bacteria, 'r', 'filled'); % Original data
fittedBacteria(:,1) = logisticFunFixedN0(fitParams, time);
plot(time, fittedBacteria(:,1), 'b-', 'LineWidth', 2); % Fitted curve
xlabel('Time');
ylabel('Bacteria Count');
title('Logistic Growth Model with Confidence Intervals');
legend('Data', 'Fitted Curve');
grid on;
hold off;
ylim([1e5 1e9]);
yticks([1e5,1e6, 1e7,1e8,1e9]);
set(gca,'FontSize',18);
set(gca,'YScale','log');


%%

bacteria = mean(CBA18,2);
% Fixed Initial Population (N0)
N0 = bacteria(1); % Set N0 to the first data point
% Logistic Growth Function with Fixed N0
logisticFunFixedN0 = @(params, t) params(1) ./ (1 + ((params(1) - N0) / N0) .* exp(-params(2) * t));
% Initial Guess for Parameters: [Carrying capacity K, Growth rate r]
initialParams = [max(bacteria), 0.1];
% Fit the model using non-linear regression with lsqcurvefit
[fitParams, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(logisticFunFixedN0, initialParams, time_no_phage, bacteria);
% Extract fitted parameters
K = fitParams(1);      % Carrying capacity
r = fitParams(2);      % Growth rate
% Calculate standard errors and confidence intervals
ci = nlparci(fitParams, residual, 'jacobian', jacobian); % 95% confidence intervals

growthrate(2) = r;
g_low(2) = ci(2,1);
g_high(2) = ci(2,2);
fittedBacteria(:,2) = logisticFunFixedN0(fitParams, time);

% Display Results
fprintf('Fitted Parameters with Confidence Intervals:\n');
fprintf('Carrying Capacity (K): %.3f (95%% CI: [%.3f, %.3f])\n', K, ci(1,1), ci(1,2));
fprintf('Growth Rate (r): %.3f (95%% CI: [%.3f, %.3f])\n', r, ci(2,1), ci(2,2));

% Plot the results
figure;
hold on;
scatter(time, bacteria, 'r', 'filled'); % Original data
fittedBacteria(:,2) = logisticFunFixedN0(fitParams, time);
plot(time, fittedBacteria(:,2), 'b-', 'LineWidth', 2); % Fitted curve
xlabel('Time');
ylabel('Bacteria Count');
title('Logistic Growth Model with Confidence Intervals');
legend('Data', 'Fitted Curve');
grid on;
hold off;
ylim([1e5 1e9]);
yticks([1e5,1e6, 1e7,1e8,1e9]);
set(gca,'FontSize',18);
set(gca,'YScale','log');
%%

bacteria = mean(CBA38,2);
% Fixed Initial Population (N0)
N0 = bacteria(1); % Set N0 to the first data point
% Logistic Growth Function with Fixed N0
logisticFunFixedN0 = @(params, t) params(1) ./ (1 + ((params(1) - N0) / N0) .* exp(-params(2) * t));
% Initial Guess for Parameters: [Carrying capacity K, Growth rate r]
initialParams = [max(bacteria), 0.1];
% Fit the model using non-linear regression with lsqcurvefit
[fitParams, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(logisticFunFixedN0, initialParams, time_no_phage, bacteria);
% Extract fitted parameters
K = fitParams(1);      % Carrying capacity
r = fitParams(2);      % Growth rate
% Calculate standard errors and confidence intervals
ci = nlparci(fitParams, residual, 'jacobian', jacobian); % 95% confidence intervals

growthrate(3) = r;
g_low(3) = ci(2,1);
g_high(3) = ci(2,2);
fittedBacteria(:,3) = logisticFunFixedN0(fitParams, time);


% Display Results
fprintf('Fitted Parameters with Confidence Intervals:\n');
fprintf('Carrying Capacity (K): %.3f (95%% CI: [%.3f, %.3f])\n', K, ci(1,1), ci(1,2));
fprintf('Growth Rate (r): %.3f (95%% CI: [%.3f, %.3f])\n', r, ci(2,1), ci(2,2));

% Plot the results
figure;
hold on;
scatter(time, bacteria, 'r', 'filled'); % Original data
fittedBacteria(:,3) = logisticFunFixedN0(fitParams, time);
plot(time, fittedBacteria(:,3), 'b-', 'LineWidth', 2); % Fitted curve
xlabel('Time');
ylabel('Bacteria Count');
title('Logistic Growth Model with Confidence Intervals');
legend('Data', 'Fitted Curve');
grid on;
hold off;
ylim([1e5 1e9]);
yticks([1e5,1e6, 1e7,1e8,1e9]);
set(gca,'FontSize',18);
set(gca,'YScale','log');
%%


bacteria = mean(H100,2);
% Fixed Initial Population (N0)
N0 = bacteria(1); % Set N0 to the first data point
% Logistic Growth Function with Fixed N0
logisticFunFixedN0 = @(params, t) params(1) ./ (1 + ((params(1) - N0) / N0) .* exp(-params(2) * t));
% Initial Guess for Parameters: [Carrying capacity K, Growth rate r]
initialParams = [max(bacteria), 0.7];
% Fit the model using non-linear regression with lsqcurvefit
[fitParams, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(logisticFunFixedN0, initialParams, time_no_phage, bacteria);
% Extract fitted parameters
K = fitParams(1);      % Carrying capacity
r = fitParams(2);      % Growth rate
% Calculate standard errors and confidence intervals
ci = nlparci(fitParams, residual, 'jacobian', jacobian); % 95% confidence intervals

growthrate(4) = r;
g_low(4) = ci(2,1);
g_high(4) = ci(2,2);
fittedBacteria(:,4) = logisticFunFixedN0(fitParams, time);



% Display Results
fprintf('Fitted Parameters with Confidence Intervals:\n');
fprintf('Carrying Capacity (K): %.3f (95%% CI: [%.3f, %.3f])\n', K, ci(1,1), ci(1,2));
fprintf('Growth Rate (r): %.3f (95%% CI: [%.3f, %.3f])\n', r, ci(2,1), ci(2,2));

% Plot the results
figure;
hold on;
scatter(time, bacteria, 'r', 'filled'); % Original data
fittedBacteria(:,4) = logisticFunFixedN0(fitParams, time);
plot(time, fittedBacteria(:,4), 'b-', 'LineWidth', 2); % Fitted curve
xlabel('Time');
ylabel('Bacteria Count');
title('Logistic Growth Model with Confidence Intervals');
legend('Data', 'Fitted Curve');
grid on;
hold off;
ylim([1e5 1e9]);
yticks([1e5,1e6, 1e7,1e8,1e9]);
set(gca,'FontSize',18);
set(gca,'YScale','log');

%%

bacteria = mean(PSA13_15,2);
% Fixed Initial Population (N0)
N0 = bacteria(1); % Set N0 to the first data point
% Logistic Growth Function with Fixed N0
logisticFunFixedN0 = @(params, t) params(1) ./ (1 + ((params(1) - N0) / N0) .* exp(-params(2) * t));
% Initial Guess for Parameters: [Carrying capacity K, Growth rate r]
initialParams = [max(bacteria), 0.7];
% Fit the model using non-linear regression with lsqcurvefit
[fitParams, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(logisticFunFixedN0, initialParams, time_no_phage, bacteria);
% Extract fitted parameters
K = fitParams(1);      % Carrying capacity
r = fitParams(2);      % Growth rate
% Calculate standard errors and confidence intervals
ci = nlparci(fitParams, residual, 'jacobian', jacobian); % 95% confidence intervals

growthrate(5) = r;
g_low(5) = ci(2,1);
g_high(5) = ci(2,2);
fittedBacteria(:,5) = logisticFunFixedN0(fitParams, time);


% Display Results
fprintf('Fitted Parameters with Confidence Intervals:\n');
fprintf('Carrying Capacity (K): %.3f (95%% CI: [%.3f, %.3f])\n', K, ci(1,1), ci(1,2));
fprintf('Growth Rate (r): %.3f (95%% CI: [%.3f, %.3f])\n', r, ci(2,1), ci(2,2));

% Plot the results
figure;
hold on;
scatter(time, bacteria, 'r', 'filled'); % Original data
fittedBacteria(:,5) = logisticFunFixedN0(fitParams, time);
plot(time, fittedBacteria(:,5), 'b-', 'LineWidth', 2); % Fitted curve
xlabel('Time');
ylabel('Bacteria Count');
title('Logistic Growth Model with Confidence Intervals');
legend('Data', 'Fitted Curve');
grid on;
hold off;
ylim([1e5 1e9]);
yticks([1e5,1e6, 1e7,1e8,1e9]);
set(gca,'FontSize',18);
set(gca,'YScale','log');

%%


load('./../../community/results/v25.mat','chain');


seivd_r = median(chain(40001:50000,28:32));
dev_seivd_r = std(chain(40001:50000,28:32));


     r_onestep = [ 0.18702      0.21341      0.21843      0.27922      0.25593];
     r_onestep_std = [    0.048525     0.066347     0.046714     0.058899     0.052958];


%% load SEIVD fits

load('./../../community/fits_seivd.mat');

%% plots



figure(2)
subplot(2,3,1)
errorbar(time,mean(host1,2),std(host1'),'LineWidth',1,'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',blue1,'MarkeredgeColor','k',Color=blue1 );hold on;
hold on;
errorbar(time_no_phage,mean(CBA4,2),std(CBA4'), 'LineWidth',1,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k',Color=blue1 )
plot(time_no_phage,fittedBacteria(:,1),'k--',LineWidth=1.5);hold on;
plot(t_after,S_after(:,1),'--',LineWidth=1,Color=blue1,Marker='none');


N0 = mean(host1(1,:)); % Initial bacterial concentration
r_dummy = seivd_r(1); % Growth rate (per unit time)
bacteria_exp = N0 * exp(r_dummy * time);
plot(time,bacteria_exp,'--',Color=[0.8500 0.3250 0.0980]);
N0 = mean(CBA4(1,:),2); % Initial bacterial concentration
r_dummy = growthrate(1); % Growth rate (per unit time)
bacteria_exp = N0 * exp(r_dummy * time_no_phage);
plot(time_no_phage,bacteria_exp,'--',Color=[0.4660 0.6740 0.1880]);

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
errorbar(time_no_phage,mean(CBA18,2),std(CBA18'), 'LineWidth',1,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k',Color=blue2 )
plot(time_no_phage,fittedBacteria(:,2),'k--',LineWidth=1.5);hold on;
plot(t_after,S_after(:,2),'--',LineWidth=1,Color=blue2,Marker='none');

N0 = mean(host2(1,:)); % Initial bacterial concentration
r_dummy = seivd_r(2); % Growth rate (per unit time)
bacteria_exp = N0 * exp(r_dummy * time);
plot(time,bacteria_exp,'--',Color=[0.8500 0.3250 0.0980]);
N0 = mean(CBA18(1,:),2); % Initial bacterial concentration
r_dummy = growthrate(2); % Growth rate (per unit time)
bacteria_exp = N0 * exp(r_dummy * time_no_phage);
plot(time_no_phage,bacteria_exp,'--',Color=[0.4660 0.6740 0.1880]);

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
errorbar(time_no_phage,mean(CBA38,2),std(CBA38'), 'LineWidth',1,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k',Color=blue3 )
plot(time_no_phage,fittedBacteria(:,3),'k--',LineWidth=1.5);hold on;
plot(t_after,S_after(:,3),'--',LineWidth=1,Color=blue3,Marker='none');

N0 = mean(host3(1,:)); % Initial bacterial concentration
r_dummy = seivd_r(3); % Growth rate (per unit time)
bacteria_exp = N0 * exp(r_dummy * time);
plot(time,bacteria_exp,'--',Color=[0.8500 0.3250 0.0980]);
N0 = mean(CBA38(1,:),2); % Initial bacterial concentration
r_dummy = growthrate(3); % Growth rate (per unit time)
bacteria_exp = N0 * exp(r_dummy * time_no_phage);
plot(time_no_phage,bacteria_exp,'--',Color=[0.4660 0.6740 0.1880]);


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
errorbar(time_no_phage,mean(H100,2),std(H100'), 'LineWidth',1,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k',Color=blue4)
plot(time_no_phage,fittedBacteria(:,4),'k--',LineWidth=1.5);hold on;
plot(t_after,S_after(:,4),'--',LineWidth=1,Color=blue4,Marker='none');

N0 = mean(host4(1,:)); % Initial bacterial concentration
r_dummy = seivd_r(4); % Growth rate (per unit time)
bacteria_exp = N0 * exp(r_dummy * time);
plot(time,bacteria_exp,'--',Color=[0.8500 0.3250 0.0980]);


N0 = mean(H100(1,:),2); % Initial bacterial concentration
r_dummy = growthrate(4); % Growth rate (per unit time)
bacteria_exp = N0 * exp(r_dummy * time_no_phage);
plot(time_no_phage,bacteria_exp,'--',Color=[0.4660 0.6740 0.1880]);


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
errorbar(time_no_phage,mean(PSA13_15,2),std(PSA13_15'), 'LineWidth',1,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k',Color=blue5)
plot(time_no_phage, fittedBacteria(:,5),'k--',LineWidth=1.5);hold on;
plot(t_after,S_after(:,5),'--',LineWidth=1,Color=blue5,Marker='none');


N0 = mean(host5(1,:)); % Initial bacterial concentration
r_dummy = seivd_r(5); % Growth rate (per unit time)
bacteria_exp = N0 * exp(r_dummy * time);
plot(time,bacteria_exp,'--',Color=[0.8500 0.3250 0.0980]);

N0 = mean(PSA13_15(1,:),2); % Initial bacterial concentration
r_dummy = growthrate(5); % Growth rate (per unit time)
bacteria_exp = N0 * exp(r_dummy * time_no_phage);
plot(time_no_phage,bacteria_exp,'--',Color=[0.4660 0.6740 0.1880]);

ylim([1e5 1e9]);
yticks([1e5,1e6, 1e7,1e8,1e9]);set(gca,'FontSize',18);
set(gca,'YScale','log');
xlabel('Time (hrs)');
xticks(0:2:16);
title('PSA 13-15');
set(gca,'fontname','times');
xlim([0 16]);
xticks(0:2:16);



% subplot(2,3,6)
% plot(time,mean(host1+host2+host3+host4+host5,2),'LineWidth',2,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor',[1,0,0],Color=[1,0,0]);hold on;
% hold on;
% plot(time,mean(CBA4+CBA18+CBA38+ H100+PSA13_15,2), 'LineWidth',2,'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',[0,0,1],Color=[0,0,1] )
% set(gca,'FontSize',18);
% set(gca,'YScale','log');
% xlabel('Time (hrs)');
% xticks(0:2:16);
% title('All hosts');
% set(gca,'fontname','times');


subplot(2,3,6)

hold on;
errorbar((1:5)-0.1,growthrate,growthrate-g_low,growthrate-g_low , 'Marker','square',MarkerSize = 8, MarkerFaceColor='w',Color='k',LineStyle='none',LineWidth=2);
errorbar((1:5),seivd_r,dev_seivd_r,'Marker','o',MarkerSize = 8, MarkerFaceColor='b',MarkerEdgeColor='k'    ,Color='b',LineStyle='none',LineWidth=2);
xticks(1:5);
errorbar((1:5)+0.1,[0.19, 0.245,0.22,0.28,0.25],[0.019,0.021,0.030,0.052,0.044],'Marker','^',MarkerEdgeColor='k'  ,MarkerSize = 8,MarkerFaceColor=[255, 219, 88]./255, MarkerEdgeColor = [255, 219, 88]./255, Color =[255, 219, 88]./255,  LineStyle='none',LineWidth=2)
%errorbar((1:5)+0.15, r_onestep,r_onestep_std, 'o','MarkerSize',8,'MarkerFaceColor',color_green,'Color',color_green,'LineWidth',2)
xticklabels({'CBA4','CBA18','CBA38','PSA H100','PSA 13-15'})
set(gca,'FontSize',18);
legend('No phage control (community): inferred','Phage infected hosts (community): inferred','Single host (doubling time experiments)');
ylabel({'Bacterial growth rates' ; '(cells/hr)'})
xlim([0.5,5.5]);
ylim([-0.2, 0.7]);
yticks(0:0.1:0.7);
grid on;
set(gca,'fontname','times');
box on;
