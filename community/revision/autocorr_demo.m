clear all;
clc;

skip  = 10 ; % from chain thinning.
% if not in use set to 1.
limit = 0.2;




load("v25.mat",'chain');
chain_2 = chain(20000:skip:end,:);

fs = 22;


%% for beta

fig1 = figure;


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

han=axes(fig1,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'\beta (virions/cell)');
title(han,'For \beta (virions/cell)');
set(gca,'FontSize',fs)


%% for phi 


fig2 = figure;


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

han=axes(fig2,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'\beta (virions/cell)');
set(gca,'FontSize',fs)


%% for tau 


fig3 = figure;


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

han=axes(fig3,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'\beta (virions/cell)');
set(gca,'FontSize',fs)


%% for r and Dc

fig4 = figure;


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





han=axes(fig4,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'\beta (virions/cell)');
set(gca,'FontSize',fs)


