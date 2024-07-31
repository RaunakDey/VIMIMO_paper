

%% load the dataset

color = [111,193,157]./255;

color1 = [76,132,147]./255;
color2 = [217,76,33]./255;


chain = chain(5000:end,:);
burn = 1;

chain_1 = chain(1:end/2,1:5);
chain_2 = chain(end/2+1:end,1:5);
mean_chain_1 = mean(chain_1);
mean_chain_2 = mean(chain_2);
N = length(chain_1)/2;
M = 2;

grand_mean = (mean_chain_1 + mean_chain_2 )/2;
B = N/(M-1) * ((mean_chain_1 - grand_mean).^2 + (mean_chain_2 - grand_mean).^2) ;
v_1 = (std(chain_1)).^2 ;
v_2 = (std(chain_2)).^2 ;
W = (1/M) * (v_1 + v_2) ;
R =  ( (N-1)/N * W + (B/N) ) ./W;

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

subplot(5,2,1)
plot(log(chain(1:end/2,2))./log(10), Color=color1);  hold on;  plot(log(chain(end/2+1:end,2))./log(10),Color= color2); 
ylabel('log \phi (mL/hr)');
set(gca, 'FontSize', 14,'FontName','Times');
xlim([0 length(chain_1)]);
title('R_{geweke} = '+string(R_geweke(2)) +  ';   R_{GL} = '+string(R(2)) );

subplot(5,2,2)
h = histogram(log(chain(1:end/2,2))./log(10),"NumBins",30,'Normalization','pdf');
bins = 0.5*(h.BinEdges(2:end) +  h.BinEdges(1:end-1));
values = h.Values;
plot(bins,values,'LineWidth',2,'Color',color)


xlabel('log \phi (mL/hr)');ylabel('PDF') ; 
set(gca, 'FontSize', 14,'FontName','Times');
hold on;
chain(:,2) = log(chain(:,2))/log(10) ;
xaxis = linspace(log( params{2,1}{1,3})/log(10), log(params{2,1}{1,4})/log(10), 200   );
% I put 4 as I remember the log-scale conversion, recheck later.
plot(xaxis, gaussian( xaxis,  log( results.prior(2,1)) /log(10),0.4 ), 'k',LineWidth=2);
xlim([xaxis(1) xaxis(end)]);




subplot(5,2,3)
plot(chain(1:end/2,3), Color= color1);  hold on;  plot(chain(end/2+1:end,3),Color=color2); 
ylabel('\tau (hr)');
set(gca, 'FontSize', 14,'FontName','Times');
xlim([0 length(chain_1)]);
title('R_{geweke} = '+string(R_geweke(3)) +  ';   R_{GL} = '+string(R(3)) );

subplot(5,2,4)
h = histogram(chain(:,3),"NumBins",30,'Normalization','pdf');
bins = 0.5*(h.BinEdges(2:end) +  h.BinEdges(1:end-1));
values = h.Values;
plot(bins,values,'LineWidth',2,'Color',color)

xlabel('\tau (hr)') ; ylabel('PDF');
set(gca, 'FontSize', 14,'FontName','Times');
hold on;

xaxis = linspace( params{3,1}{1,3}, params{3,1}{1,4}, 200   );
plot(xaxis, gaussian( xaxis,   results.prior(3,1),results.prior(3,2)),'k',LineWidth=2);
xlim([xaxis(1) xaxis(end)]);




subplot(5,2,5)
plot(chain(1:end/2,4),Color= color1);  hold on;  plot(chain(end/2+1:end,4),Color=color2); 
ylabel('\beta (virions/cell)');
set(gca, 'FontSize', 14,'FontName','Times');
xlim([0 length(chain_1)]);
title('R_{geweke} = '+string(R_geweke(4)) +  ';   R_{GL} = '+string(R(4)) );

subplot(5,2,6)
h = histogram(chain(:,4),"NumBins",30,'Normalization','pdf');
bins = 0.5*(h.BinEdges(2:end) +  h.BinEdges(1:end-1));
values = h.Values;
plot(bins,values,'LineWidth',2,'Color',color)


xlabel('\beta (virions/cell)') ; ylabel('PDF');
set(gca, 'FontSize', 14,'FontName','Times');
hold on;
xaxis = linspace( params{4,1}{1,3}, params{4,1}{1,4}, 200   );
plot(xaxis, gaussian( xaxis,   results.prior(4,1),results.prior(4,2)),'k',LineWidth=2);

subplot(5,2,7)
plot(chain(1:end/2,5), Color=color1);  hold on;  plot(chain(end/2+1:end,5),Color=color2);
xlabel('MCMC step');ylabel('N_E');
set(gca, 'FontSize', 14,'FontName','Times');
xlim([0 length(chain_1)]);
title('R_{geweke} = '+string(R_geweke(5)) +  ';   R_{GL} = '+string(R(5)) );




subplot(5,2,8)
h = histogram(chain(:,5),"NumBins",30,'Normalization','pdf');
bins = 0.5*(h.BinEdges(2:end) +  h.BinEdges(1:end-1));
values = h.Values;
plot(bins,values,'LineWidth',2,'Color',color)
xlim([xaxis(1) xaxis(end)]);

xlabel('N_E') ; ylabel('PDF');
set(gca, 'FontSize', 14,'FontName','Times');
hold on;
xaxis = linspace( params{5,1}{1,3}, params{5,1}{1,4}, 200   );
plot(xaxis, gaussian( xaxis,   results.prior(5,1),results.prior(5,2)),'k',LineWidth=2);
xlim([xaxis(1) xaxis(end)]);





subplot(5,2,9)
plot(chain(1:end/2,1), Color=color1);  hold on;  plot(chain(end/2+1:end,1),Color=color2);
xlabel('MCMC step');ylabel('r (cells/hr)');
set(gca, 'FontSize', 14,'FontName','Times');
xlim([0 length(chain_1)]);
title('R_{geweke} = '+string(R_geweke(1)) +  ';   R_{GL} = '+string(R(1)) );




subplot(5,2,10)
h = histogram(chain(:,1),"NumBins",30,'Normalization','pdf');
bins = 0.5*(h.BinEdges(2:end) +  h.BinEdges(1:end-1));
values = h.Values;
plot(bins,values,'LineWidth',2,'Color',color)
xlim([xaxis(1) xaxis(end)]);

xlabel('r (cells/hr)') ; ylabel('PDF');
set(gca, 'FontSize', 14,'FontName','Times');
hold on;
xaxis = linspace( params{1,1}{1,3}, params{1,1}{1,4}, 200   );
plot(xaxis, gaussian( xaxis,   results.prior(1,1),results.prior(1,2)),'k',LineWidth=2);
xlim([xaxis(1) xaxis(end)]);



clc;
clear all;