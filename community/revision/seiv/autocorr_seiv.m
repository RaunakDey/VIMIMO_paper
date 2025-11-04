clear all;
clc;



load("cba_seiv_2.mat",'chain_final');
chain_cba = chain_final; 
load("psa_seiv_2.mat",'chain_final');
chain_psa = chain_final;

chain_joined(:,1:5) = chain_cba(:,1:5);
chain_joined(:,6:9) = chain_psa(:,1:4);
chain_joined(:,10:14) = chain_cba(:,6:10);
chain_joined(:,15:18) = chain_psa(:,5:8);
chain_joined(:,19:23) = chain_cba(:,11:15);
chain_joined(:,24:27) = chain_psa(:,9:12);
chain_joined(:,28:30) = chain_cba(:,16:18);
chain_joined(:,31:32) = chain_psa(:,13:14);



skip  = 10 ; % from chain thinning.
% if not in use set to 1.
limit = 0.2;





chain_2 = chain_joined(20000:skip:end,:);

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

if i==5
ylabel('Autocorrelation');
end
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

if i==5
ylabel('Autocorrelation');
end

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

if i == 5
ylabel('Autocorrelation');
end

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
    subplot(5,1,i)
[acf,lags] = autocorr(chain_2(:,i+27),NumLags=200);
bh = bar(acf);
bh.FaceColor = [62,137,168]./255;


yline(limit,'--',Color='k',LineWidth=2)
yline(-limit,'--',Color='k',LineWidth=2)

if i == 3
ylabel('Autocorrelation')
end

make_title(i);
set(gca,'FontSize',fs)
end








han=axes(fig4,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'\beta (virions/cell)');
set(gca,'FontSize',fs)


