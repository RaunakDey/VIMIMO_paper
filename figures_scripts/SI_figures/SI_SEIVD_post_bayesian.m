clc;
clear all;


%% load data

load('./../../community/results/v25.mat');

%% priors and posteriors
chain2 = chain;


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
smoothHistogram(chain(transient_id:end,i),6,'b');
%smoothHistogram(chain2(transient_id:end,i),6,'r');
%histogram(chain(transient_id:end,i),"Normalization","pdf",DisplayStyle="stairs");
%histogram(chain2(transient_id:end,i), "Normalization","pdf",DisplayStyle="stairs");

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
%smoothHistogram(chain2(transient_id:end,i),6,'r');

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
%smoothHistogram(chain2(transient_id:end,i),6,'r');

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
%smoothHistogram(chain2(transient_id:end,i),6,'r');

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
%smoothHistogram(chain2(transient_id:end,i),6,'r');

xlabel('Dc');
ylabel('PDF');
xlim([min(chain(transient_id:end,i)) max(chain(transient_id:end,i)) ]);
set(gca,'FontSize',20);
end

%% trace plots



fig1 = figure;
load('./../../community/results/v23.mat');
chain_1 = chain(20000:10:50000,:);

load('./../../community/results/v25.mat');
chain_2 = chain(20000:10:50000,:);



transient_id = 1; %already took out the burnout.
 

color1 = [76,132,147]./255;
color2 = [217,76,33]./255;

for i = 1:9
subplot(5,2,i)
plot(chain_1(transient_id:end,i),'Color',color1 );
xlim([0 3000])
hold on;
end


for i = 1:9
subplot(5,2,i)
plot(chain_2(transient_id:end,i),'Color',color2);
ylim([  min(min(chain_1(transient_id:end,i)),min(chain_2(transient_id:end,i)))-0.1,  max(max(chain_1(transient_id:end,i)),max(chain_2(transient_id:end,i)))+0.1 ]);
xlim([0 3000])
make_title(i);
end

han=axes(fig1,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'\beta (virions/cell)');
xlabel(han,'steps');
%title(han,'yourTitle');
set(gca,'FontSize',18)






fig2 = figure;
for i = 1:9
subplot(5,2,i)
plot(chain_1(transient_id:end,9+i), 'Color',color1 );
xlim([0 3000])
hold on;
end

for i = 1:9
subplot(5,2,i)
plot(chain_2(transient_id:end,9+i),'Color',color2);
ylim([  min(min(chain_1(transient_id:end,9+i)),min(chain_2(transient_id:end,9+i)))-0.1,  max(max(chain_1(transient_id:end,9+i)),max(chain_2(transient_id:end,9+i)))+0.1 ]);
xlim([0 3000])
make_title(i);
end

han=axes(fig2,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'log_{10} \phi (ml/hr) ');
xlabel(han,'steps');
%title(han,'yourTitle');
set(gca,'FontSize',18)



fig3 = figure;
for i = 1:9

subplot(5,2,i)
plot(chain_1(transient_id:end,18+i),'Color',color1);
xlim([0 3000])
hold on;

end

for i = 1:9
subplot(5,2,i)
plot(chain_2(transient_id:end,18+i),'Color',color2);
ylim([  min(min(chain_1(transient_id:end,18+i)),min(chain_2(transient_id:end,18+i)) -0.1) ,  max(max(chain_1(transient_id:end,18+i)),max(chain_2(transient_id:end,18+i)) +0.1)]);
xlim([0 3000])
make_title(i);
end


han=axes(fig3,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'\tau (hr) ');
xlabel(han,'steps');
%title(han,'yourTitle');
set(gca,'FontSize',18)



fig4 = figure;
 
for i = 1:5
subplot(3,2,i)
plot(chain_1(transient_id:end,i+27), 'Color',color1 );
xlim([0 3000])
hold on;
end


for i = 1:5
subplot(3,2,i)
plot(chain_2(transient_id:end,i+27),'Color',color2);
ylim([  min(min(chain_1(transient_id:end,i+27)),min(chain_2(transient_id:end,i+27))),  max(max(chain_1(transient_id:end,i+27)),max(chain_2(transient_id:end,i+27))) ]);
xlim([0 3000])
make_title_host(i);
end

han=axes(fig4,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'r (cells/hr)');
xlabel(han,'steps');
%title(han,'yourTitle');
set(gca,'FontSize',18)

