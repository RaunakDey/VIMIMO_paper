chain = chain(transient_id:end,:);
transient_id = 1;
limit = 0.2;
%% covar



% figure()
% subplot(1,3,1)
% i= 2;
% j= 3;
% pearson = corrcoef(chain(transient_id:end,i),chain(transient_id:end,j));
% pearson_corr = pearson(1,2);
% 
% if pearson_corr>0
% color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
% else
% color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
% end
% 
% plot(chain(transient_id:end,j),chain(transient_id:end,i),'.','Color',color_input);
% set(gca,'FontSize',14);
% axis square;
% ylabel('\phi (mL/hr)');
% xlabel('\tau (hr)');
% 
% 
% 
% subplot(1,3,2)
% i= 3;
% j= 4;
% pearson = corrcoef(chain(transient_id:end,i),chain(transient_id:end,j));
% pearson_corr = pearson(1,2);
% 
% if pearson_corr>0
% color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
% else
% color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
% end
% 
% plot(chain(transient_id:end,j),chain(transient_id:end,i),'.','Color',color_input);
% set(gca,'FontSize',14);
% axis square;
% 
% ylabel('\tau (hr)');
% xlabel('\beta (virions/cells)');
% 
% 
% subplot(1,3,3)
% i= 4;
% j= 2;
% pearson = corrcoef(chain(transient_id:end,i),chain(transient_id:end,j));
% pearson_corr = pearson(1,2);
% 
% if pearson_corr>0
% color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
% else
% color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
% end
% 
% plot(chain(transient_id:end,j),chain(transient_id:end,i),'.','Color',color_input);
% set(gca,'FontSize',14);
% axis square;
% ylabel('\tau (hr)');
% xlabel('\phi (mL/hr)');

%% load the dataset

color = [111,193,157]./255;

color1 = [76,132,147]./255;
color2 = [217,76,33]./255;


%chain = chain(5000:end,:);
%chain = chain(1:end,:);

burn = 1;

chain_1 = chain(1:end/2,1:5);
chain_2 = chain(end/2+1:end,1:5);
mean_chain_1 = mean(chain_1);
mean_chain_2 = mean(chain_2);
N = length(chain_1);
M = 2;

grand_mean = (mean_chain_1 + mean_chain_2 )/2;
B = N/(M-1) * ((mean_chain_1 - grand_mean).^2 + (mean_chain_2 - grand_mean).^2) ;
v_1 = (std(chain_1)).^2 ;
v_2 = (std(chain_2)).^2 ;
W = (1/M) * (v_1 + v_2) ;
R =  sqrt(( (N-1)/N * W + (B/N) ) ./W);

% for i = 1:5
%     if R(i) <1
%         R(i) = 1.0;
%     end
% end

store = chainstats(chain(1:end,:));
R_geweke = store(:,5);

clear min;
clear max;


figure

subplot(4,3,1)
plot(log(chain(1:end/2,2))./log(10), Color=color1);  hold on;  plot(log(chain(end/2+1:end,2))./log(10),Color= color2); 
ylabel('log \phi (mL/hr)');
set(gca, 'FontSize', 14,'FontName','Times');
xlim([0 length(chain_1)]);
%title('R_{geweke} = '+string(R_geweke(2)) +  ';   R_{GL} = '+string(R(2)) );
title('R_{GL} = '+string(R(2)) );


subplot(4,3,2)
% h = histogram(log(chain(1:end/2,2))./log(10),"NumBins",30,'Normalization','pdf');
% bins = 0.5*(h.BinEdges(2:end) +  h.BinEdges(1:end-1));
% values = h.Values;
% plot(bins,values,'LineWidth',2,'Color',color)

[fp,xfp] = kde(log(chain(1:end/2,2))./log(10));
plot(xfp,fp,LineWidth=2,Color='k');
clear xfp fp;



xlabel('log \phi (mL/hr)');ylabel('PDF') ; 
set(gca, 'FontSize', 14,'FontName','Times');
hold on;
chain(:,2) = log(chain(:,2))/log(10) ;
xaxis = linspace(log( params{2,1}{1,3})/log(10), log(params{2,1}{1,4})/log(10), 200   );
% I put 4 as I remember the log-scale conversion, recheck later.
plot(xaxis, gaussian( xaxis,  log( results.prior(2,1)) /log(10),0.4 ), 'k',LineWidth=2,LineStyle='--');
xlim([xaxis(1) xaxis(end)]);
xlim([-10 -5]);




subplot(4,3,3)
[acf,lags] = autocorr(chain(:,2),NumLags=200);
bh = bar(acf);
bh.FaceColor = [0 0 0];
ylim([-0.22 1.01]);
xlim([0 200]);
xticks([0,100,200]);
yline(limit,'--',Color='k',LineWidth=2)
yline(-limit,'--',Color='k',LineWidth=2)
ylabel('Autocorrelation');

%%

subplot(4,3,4)
plot(chain(1:end/2,3), Color= color1);  hold on;  plot(chain(end/2+1:end,3),Color=color2); 
ylabel('\tau (hr)');
set(gca, 'FontSize', 14,'FontName','Times');
xlim([0 length(chain_1)]);
%title('R_{geweke} = '+string(R_geweke(3)) +  ';   R_{GL} = '+string(R(3)) );
title('R_{GL} = '+string(R(3)) );

subplot(4,3,5)
% h = histogram(chain(:,3),"NumBins",30,'Normalization','pdf');
% bins = 0.5*(h.BinEdges(2:end) +  h.BinEdges(1:end-1));
% values = h.Values;
% plot(bins,values,'LineWidth',2,'Color',color)

[fp,xfp] = kde(chain(:,3));
plot(xfp,fp,LineWidth=2,Color='k');
clear xfp fp;


xlabel('\tau (hr)') ; ylabel('PDF');
set(gca, 'FontSize', 14,'FontName','Times');
hold on;

xaxis = linspace( params{3,1}{1,3}, params{3,1}{1,4}, 200   );
plot(xaxis, gaussian( xaxis,   results.prior(3,1),results.prior(3,2)),'k',LineWidth=2,LineStyle='--');
xlim([xaxis(1) xaxis(end)]);
xlim([0.5 5]);


subplot(4,3,6)
[acf,lags] = autocorr(chain(:,3),NumLags=200);
bh = bar(acf);
bh.FaceColor = [0 0 0];
ylim([-0.22 1.01]);
xlim([0 200]);
xticks([0,100,200]);
yline(limit,'--',Color='k',LineWidth=2)
yline(-limit,'--',Color='k',LineWidth=2)
ylabel('Autocorrelation');


%%

subplot(4,3,7)
plot(chain(1:end/2,4),Color= color1);  hold on;  plot(chain(end/2+1:end,4),Color=color2); 
ylabel('\beta (virions/cell)');
set(gca, 'FontSize', 14,'FontName','Times');
xlim([0 length(chain_1)]);
%title('R_{geweke} = '+string(R_geweke(4)) +  ';   R_{GL} = '+string(R(4)) );
title('R_{GL} = '+string(R(4)) );

subplot(4,3,8)
% h = histogram(chain(:,4),"NumBins",30,'Normalization','pdf');
% bins = 0.5*(h.BinEdges(2:end) +  h.BinEdges(1:end-1));
% values = h.Values;
% plot(bins,values,'LineWidth',2,'Color',color)

[fp,xfp] = kde(chain(:,4));
plot(xfp,fp,LineWidth=2,Color='k');
clear xfp fp;

xlabel('\beta (virions/cell)') ; ylabel('PDF');
set(gca, 'FontSize', 14,'FontName','Times');
hold on;
xaxis = linspace( params{4,1}{1,3}, params{4,1}{1,4}, 200   );
plot(xaxis, gaussian( xaxis,   results.prior(4,1),results.prior(4,2)),'k',LineWidth=2,LineStyle='--');
xlim([0.5 600]);


subplot(4,3,9)
[acf,lags] = autocorr(chain(:,4),NumLags=200);
bh = bar(acf);
bh.FaceColor = [0 0 0];
ylim([-0.22 1.01]);
xlim([0 200]);
xticks([0,100,200]);
yline(limit,'--',Color='k',LineWidth=2)
yline(-limit,'--',Color='k',LineWidth=2)
ylabel('Autocorrelation');



%%




subplot(4,3,10)
plot(chain(1:end/2,5), Color=color1);  hold on;  plot(chain(end/2+1:end,5),Color=color2);
xlabel('MCMC step');ylabel('N_E');
set(gca, 'FontSize', 14,'FontName','Times');
xlim([0 length(chain_1)]);
%title('R_{geweke} = '+string(R_geweke(5)) +  ';   R_{GL} = '+string(R(5)) );
title('R_{GL} = '+string(R(5)) );




subplot(4,3,11)
h = histogram(chain(:,5),"NumBins",30,'Normalization','pdf');
bins = 0.5*(h.BinEdges(2:end) +  h.BinEdges(1:end-1));
values = h.Values;
plot(bins,values,'LineWidth',2,'Color',color)

[fp,xfp] = kde(chain(:,5));
plot(xfp,fp,LineWidth=2,Color='k');
clear xfp fp;


xlim([xaxis(1) xaxis(end)]);

xlabel('N_E') ; ylabel('PDF');
set(gca, 'FontSize', 14,'FontName','Times');
hold on;
xaxis = linspace( params{5,1}{1,3}, params{5,1}{1,4}, 200   );
plot(xaxis, gaussian( xaxis,   results.prior(5,1),results.prior(5,2)),'k',LineWidth=2,LineStyle='--');
xlim([xaxis(1) xaxis(end)]);
xlim([5 200]);



subplot(4,3,12)
[acf,lags] = autocorr(chain(:,5),NumLags=200);
bh = bar(acf);
bh.FaceColor = [0 0 0];
ylim([-0.22 1.01]);
xlim([0 200]);
xticks([0,100,200]);
yline(limit,'--',Color='k',LineWidth=2)
yline(-limit,'--',Color='k',LineWidth=2)
ylabel('Autocorrelation');
xlabel('Lags');



%%

clc;
clear all;