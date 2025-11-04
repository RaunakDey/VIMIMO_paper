clear all;
clc;

skip  = 5 ; % from chain thinning.
% if not in use set to 1.

load("cba_seiv_2.mat",'chain_final','theta_initial','theta_std');
chain_cba = chain_final; 
theta_initial_cba = theta_initial;
theta_std_cba = theta_std;

load("psa_seiv_2.mat",'chain_final','theta_initial','theta_std');
chain_psa = chain_final;
theta_initial_psa = theta_initial;
theta_std_psa = theta_std;



chain_joined(:,1:5) = chain_cba(1:end,1:5);
chain_joined(:,6:9) = chain_psa(1:end,1:4);
chain_joined(:,10:14) = chain_cba(1:end,6:10);
chain_joined(:,15:18) = chain_psa(1:end,5:8);
chain_joined(:,19:23) = chain_cba(1:end,11:15);
chain_joined(:,24:27) = chain_psa(1:end,9:12);
chain_joined(:,28:30) = chain_cba(1:end,16:18);
chain_joined(:,31:32) = chain_psa(1:end,13:14);

mu_prior(1:5) = theta_initial_cba(1:5);
mu_prior(6:9) = theta_initial_psa(1:4);
mu_prior(10:14) = theta_initial_cba(6:10);
mu_prior(15:18) = theta_initial_psa(5:8);
mu_prior(19:23) = theta_initial_cba(11:15);
mu_prior(24:27) = theta_initial_psa(9:12);
mu_prior(28:30) = theta_initial_cba(16:18);
mu_prior(31:32) = theta_initial_psa(13:14);


sigma_prior(1:5) = theta_std_cba(1:5);
sigma_prior(6:9) = theta_std_psa(1:4);
sigma_prior(10:14) = theta_std_cba(6:10);
sigma_prior(15:18) = theta_std_psa(5:8);
sigma_prior(19:23) = theta_std_cba(11:15);
sigma_prior(24:27) = theta_std_psa(9:12);
sigma_prior(28:30) = theta_std_cba(16:18);
sigma_prior(31:32) = theta_std_psa(13:14);



fs = 22;


%% add paths
addpath(genpath('./../../'));


%% construct the prior

% all of them are gaussian priors in log or linear scales


transient_id = 45000;
% comment away if not using scaling based parameter reparameterization
chain = chain_joined;

figure(1)

for i=1:9
subplot(2,5,i)
xaxis = linspace( mu_prior(i) - 3*sigma_prior(i),  mu_prior(i) + 3*sigma_prior(i) ,100);
plot(xaxis, ...
    gaussian(xaxis,mu_prior(i),sigma_prior(i)), ...
    "Color",'k',LineWidth=2);
hold on;
smoothHistogram(chain(transient_id:end,i),20,'b');
xlim([min(chain(transient_id:end,i)) max(chain(transient_id:end,i)) ]);
xlabel('\beta (virions/cell)');
ylabel('PDF');
set(gca, 'FontSize', 15);
axis('square');

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
smoothHistogram(chain(transient_id:end,i),10,'b');
xlabel('\phi (ml/hr/cell)');
ylabel('PDF');%xlim([min(chain(transient_id:end,i)) max(chain(transient_id:end,i)) ]);
set(gca, 'FontSize', 15);
axis('square');

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
smoothHistogram(chain(transient_id:end,i),10,'b');
xlabel('\tau (hr)');
ylabel('PDF');
xlim([min(chain(transient_id:end,i)) max(chain(transient_id:end,i)) ]);
set(gca, 'FontSize', 15);
axis('square');

end




figure(4)

for i=1:5
subplot(1,5,i)
i=i+27;
xaxis = linspace( mu_prior(i) - 3*sigma_prior(i),  mu_prior(i) + 3*sigma_prior(i) ,100);
plot(xaxis, ...
    gaussian(xaxis,mu_prior(i),sigma_prior(i)), ...
    "Color",'k',LineWidth=2);
hold on;
smoothHistogram(chain(transient_id:end,i),10,'b');
xlabel('r (/hr)');
ylabel('PDF');
xlim([min(chain(transient_id:end,i)) max(chain(transient_id:end,i)) ]);
set(gca, 'FontSize', 15);
axis('square');

end

