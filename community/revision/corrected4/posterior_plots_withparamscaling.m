close all;
clear all;
clc;

%% import datasets as matlab file

%load('./from_os_prior_20000.mat');

%load('./v14-4.mat')     

load('./v27.mat');
chain2 = chain;

load('./v25.mat');





%% add paths
addpath(genpath('./../../'));


%% construct the prior

% all of them are gaussian priors in log or linear scales
%mu_prior = mcmcresults.prior(:,1);
%sigma_prior = mcmcresults.prior(:,2);

mu_prior = theta_start;
sigma_prior = theta_std;

transient_id = 30000;
%chain = repmat(theta_start,mcmcoptions.nsimu,1) + chain.*repmat(theta_std,mcmcoptions.nsimu,1);


figure(1)

for i=1:9
subplot(2,5,i)
xaxis = linspace( mu_prior(i) - 3*sigma_prior(i),  mu_prior(i) + 3*sigma_prior(i) ,100);
plot(xaxis, ...
    gaussian(xaxis,mu_prior(i),sigma_prior(i)), ...
    "Color",'k',LineWidth=2);
hold on;
%smoothHistogram(chain(transient_id:end,i),6,'b');
%smoothHistogram(chain2(transient_id:end,i),6,'r');
histogram(chain(transient_id:end,i),"Normalization","pdf",DisplayStyle="stairs");
histogram(chain2(transient_id:end,i), "Normalization","pdf",DisplayStyle="stairs");

xlim([min(chain(transient_id:end,i))-50 max(chain(transient_id:end,i))+50 ]);
xlabel('\beta (virions/cell)');
ylabel('PDF');
set(gca,'FontSize',20);
end



figure(2)
for i=1:9
subplot(2,5,i)
i=i+9;
xaxis = linspace( mu_prior(i) - 3*sigma_prior(i),  mu_prior(i) + 3*sigma_prior(i) ,100);
plot(xaxis, ...
    gaussian(xaxis,mu_prior(i),sigma_prior(i)), ...
    "Color",'k',LineWidth=2);
hold on;
smoothHistogram(chain(transient_id:end,i),6,'b');
smoothHistogram(chain2(transient_id:end,i),6,'r');

xlabel('\phi (ml/hr/cell)');
ylabel('PDF');
xlim([min(chain(transient_id:end,i))-1 max(chain(transient_id:end,i))+1 ]);
set(gca,'FontSize',20);
end




figure(3)

for i=1:9
subplot(2,5,i)
i=i+18;
xaxis = linspace( mu_prior(i) - 3*sigma_prior(i),  mu_prior(i) + 3*sigma_prior(i) ,100);
plot(xaxis, ...
    gaussian(xaxis,mu_prior(i),sigma_prior(i)), ...
    "Color",'k',LineWidth=2);
hold on;
smoothHistogram(chain(transient_id:end,i),6,'b');
smoothHistogram(chain2(transient_id:end,i),6,'r');

xlabel('\tau (hr)');
ylabel('PDF');
xlim([min(chain(transient_id:end,i))-0.5 max(chain(transient_id:end,i))+0.5 ]);
set(gca,'FontSize',20);
end




figure(4)

for i=1:5
subplot(2,5,i)
i=i+27;
xaxis = linspace( mu_prior(i) - 3*sigma_prior(i),  mu_prior(i) + 3*sigma_prior(i) ,100);
plot(xaxis, ...
    gaussian(xaxis,mu_prior(i),sigma_prior(i)), ...
    "Color",'k',LineWidth=2);
hold on;
smoothHistogram(chain(transient_id:end,i),6,'b');
smoothHistogram(chain2(transient_id:end,i),6,'r');

xlabel('r (/hr)');
ylabel('PDF');
xlim([min(chain(transient_id:end,i)) max(chain(transient_id:end,i)) ]);
set(gca,'FontSize',20);
end

for i=1:5
subplot(2,5,i+5)
i=i+32;
xaxis = linspace( mu_prior(i) - 3*sigma_prior(i),  mu_prior(i) + 3*sigma_prior(i) ,100);
plot(xaxis, ...
    gaussian(xaxis,mu_prior(i),sigma_prior(i)), ...
    "Color",'k',LineWidth=2);
hold on;
smoothHistogram(chain(transient_id:end,i),6,'b');
smoothHistogram(chain2(transient_id:end,i),6,'r');

xlabel('Dc');
ylabel('PDF');
xlim([min(chain(transient_id:end,i)) max(chain(transient_id:end,i)) ]);
set(gca,'FontSize',20);
end

