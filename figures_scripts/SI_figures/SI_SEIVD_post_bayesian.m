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

xlabel('log(\phi) (ml/hr/cell)');
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

xlabel('log(Dc) (dead cells/ml)');
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



fig5 = figure;
 
for i = 1:5
subplot(3,2,i)
plot(chain_1(transient_id:end,i+32), 'Color',color1 );
xlim([0 3000])
hold on;
end


for i = 1:5
subplot(3,2,i)
plot(chain_2(transient_id:end,i+32),'Color',color2);
ylim([  min(min(chain_1(transient_id:end,i+32)),min(chain_2(transient_id:end,i+32))),  max(max(chain_1(transient_id:end,i+32)),max(chain_2(transient_id:end,i+32))) ]);
xlim([0 3000])
make_title_host(i);
end

han=axes(fig5,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'log_{10}D_c (dead cells/mL)');
xlabel(han,'steps');
%title(han,'yourTitle');
set(gca,'FontSize',18)


%% autocorr


skip  = 10 ; % from chain thinning.
% if not in use set to 1.
limit = 0.2;




load("v25.mat",'chain');
chain_2 = chain(20000:skip:end,:);

fs = 22;



fig6 = figure;


for i = 1:9
    subplot(5,2,i)
[acf,lags] = autocorr(chain_2(:,i),NumLags=200);
bh = bar(acf);
bh.FaceColor = [211,23,24]./255;
yline(limit,'--',Color='k',LineWidth=2)
yline(-limit,'--',Color='k',LineWidth=2)

ylabel('Autocorrelation')
make_title(i);
set(gca,'FontSize',fs)
end

han=axes(fig6,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'\beta (virions/cell)');
title(han,'For \beta (virions/cell)');
set(gca,'FontSize',fs)



fig7 = figure;


for i = 1:9
    subplot(5,2,i)
[acf,lags] = autocorr(chain_2(:,i+9),NumLags=200);
bh = bar(acf);
bh.FaceColor = [15,104,82]./255;


yline(limit,'--',Color='k',LineWidth=2)
yline(-limit,'--',Color='k',LineWidth=2)

ylabel('Autocorrelation')
make_title(i);
set(gca,'FontSize',fs)
end

han=axes(fig7,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'\beta (virions/cell)');
set(gca,'FontSize',fs)


fig8 = figure;


for i = 1:9
    subplot(5,2,i)
[acf,lags] = autocorr(chain_2(:,i+18),NumLags=200);
bh = bar(acf);
bh.FaceColor = [211,119,46]./255;


yline(limit,'--',Color='k',LineWidth=2)
yline(-limit,'--',Color='k',LineWidth=2)
ylabel('Autocorrelation')
make_title(i);
set(gca,'FontSize',fs)
end

han=axes(fig8,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'\beta (virions/cell)');
set(gca,'FontSize',fs)



fig9 = figure;


for i = 1:5
    subplot(5,2,2*i-1)
[acf,lags] = autocorr(chain_2(:,i+27),NumLags=200);
bh = bar(acf);
bh.FaceColor = [62,137,168]./255;


yline(limit,'--',Color='k',LineWidth=2)
yline(-limit,'--',Color='k',LineWidth=2)
ylabel('Autocorrelation')
make_title(i);
set(gca,'FontSize',fs)
end



for i = 1:5
    subplot(5,2,2*i)
[acf,lags] = autocorr(chain_2(:,i+32),NumLags=200);
bh = bar(acf);
bh.FaceColor = [217,76,33]./255;


yline(limit,'--',Color='k',LineWidth=2)
yline(-limit,'--',Color='k',LineWidth=2)
ylabel('Autocorrelation')
make_title(i);
set(gca,'FontSize',fs)
end





han=axes(fig9,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'\beta (virions/cell)');
set(gca,'FontSize',fs)



%% convergence stats.

load('./../../community/results/v25.mat');
chain_1 = chain(10001:10:30000,:);
chain_2 = chain(30001:10:50000,:);

mean_chain_1 = mean(chain_1);
mean_chain_2 = mean(chain_2);
N=length(chain_1);
M=2;

grand_mean = (mean_chain_1 + mean_chain_2 )/2;

%chain_1 and chain_2 are the chains.

B = N/(M-1) * ((mean_chain_1 - grand_mean).^2 + (mean_chain_2 - grand_mean).^2) ;

v_1 = (std(chain_1)).^2 ;
v_2 = (std(chain_2)).^2 ;

W = (1/M) * (v_1 + v_2) ;
R_1 =  ( (N-1)/N * W + (B/N) ) ./W;



load('./../../community/results/v23.mat');
chain_1 = chain(30001:10:40000,:);
chain_2 = chain(40001:10:50000,:);

mean_chain_1 = mean(chain_1);
mean_chain_2 = mean(chain_2);
N=length(chain_1);
M=2;

grand_mean = (mean_chain_1 + mean_chain_2 )/2;

%chain_1 and chain_2 are the chains.

B = N/(M-1) * ((mean_chain_1 - grand_mean).^2 + (mean_chain_2 - grand_mean).^2) ;

v_1 = (std(chain_1)).^2 ;
v_2 = (std(chain_2)).^2 ;

W = (1/M) * (v_1 + v_2) ;
R_2 =  ( (N-1)/N * W + (B/N) ) ./W;


R =  min(R_1,R_2);

figure(100)
subplot(1,2,1)
b = bar(R_2,'FaceColor','flat');

for i = 1:9
b.CData(i,:) = [211,23,24]./255;
end

for i = 10:18
b.CData(i,:) = [15,104,82]./255;
end

for i = 19:27
b.CData(i,:) = [211,119,46]./255;
end

for i=28:32
b.CData(i,:) = [62,137,168]./255;
end

for i=33:37
b.CData(i,:) = [217,76,33]./255;
end


set(gca,'XTickLabel','');
hold on;
line([0.2, 38], [1.1, 1.1],'Color', [0.1,0.1,0.1],'LineStyle','--',LineWidth=2  );
set(gca,'FontSize',20);
xlim([0 38])
ylabel('R_{GR}')
yticks([0 0.2 0.4 0.6 0.8 0.9 1]);


set(gca, 'XTick', [5 14 23 30 35]);
set(gca, 'XTickLabel', { '\beta_{ij} ' '\phi_{ij}' '\tau_{ij}' 'r_i' '{Dc}_i'});
set(gca, 'fontname','times');


%% geeweke part
load('./../../community/results/v23.mat');
chainstats1 = chainstats(chain,mcmcresults);
geweke1 = chainstats1(:,5);

load('./../../community/results/v25.mat');
chainstats2 = chainstats(chain,mcmcresults);
geweke2 = chainstats2(:,5);


geweke_stats = max(geweke1,geweke2);

figure(100)
subplot(1,2,2)
b = bar(geweke_stats,'FaceColor','flat');

for i = 1:9
b.CData(i,:) = [211,23,24]./255;
end

for i = 10:18
b.CData(i,:) = [15,104,82]./255;
end

for i = 19:27
b.CData(i,:) = [211,119,46]./255;
end

for i=28:32
b.CData(i,:) = [62,137,168]./255;
end

for i=33:37
b.CData(i,:) = [217,76,33]./255;
end


set(gca,'XTickLabel','');
hold on;
line([0.2, 38], [0.9, 0.9],'Color', [0.1,0.1,0.1],'LineStyle','--',LineWidth=2  );
set(gca,'FontSize',20);
xlim([0 38])
ylabel('R_{geweke}')
yticks([0 0.2 0.4 0.6 0.8 0.9 1]);


set(gca, 'XTick', [5 14 23 30 35]);
set(gca, 'XTickLabel', { '\beta_{ij} ' '\phi_{ij}' '\tau_{ij}' 'r_i' '{Dc}_i'});
set(gca, 'fontname','times');


%% covariance plots.


load('./../../community/results/v23.mat');
%load('./../../community/results/v25.mat');

n = 256; % Number of colors in the colormap
blue_gradient = linspace(0,1,n)'.*[1,1,1] + (1-linspace(0,1,n))'.*[0,0,1];
red_gradient = linspace(0,1,n)'.*[1,1,1] + (1-linspace(0,1,n))'.*[1,0,0];
% Combine gradients into a colormap
blueredColormap = [red_gradient; flipud(blue_gradient)];



correlation_pearson = corrcoef(chain(30000:10:end,:));

for i = 1:37
    for j = 1:37
        if i<=j | abs(correlation_pearson(i,j))<0.5
        %if i<=j
            correlation_pearson(i,j) = 0;
        end
    end
end



% this will not be used in the paper and I will only revisit them if I have
% to.



figure
imagesc(correlation_pearson);
colorbar;
clim([-1, 1]); % Set the limits of the color axis
colormap(blueredColormap );
title("pearson correlation between posteriors");
axis("square");
set(gca,'fontname','times');
set(gca,'FontSize',24);
set(gca, 'XTick', [5 14 23 30 35]);
set(gca, 'XTickLabel', { '\beta_{ij} ' '\phi_{ij}' '\tau_{ij}' 'r_i' '{Dc}_i'});
set(gca, 'YTick', [5 14 23 30 35]);
set(gca, 'YTickLabel', { '\beta_{ij} ' '\phi_{ij}' '\tau_{ij}' 'r_i' '{Dc}_i'});



