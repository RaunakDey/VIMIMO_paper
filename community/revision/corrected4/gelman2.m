clear all;
clc;


load("v23.mat",'chain');
chain_1 = chain(20001:end,:);


load("v25.mat",'chain');
chain_2 = chain(20001:end,:);


mean_chain_1 = mean(chain_1);
mean_chain_2 = mean(chain_2);


%%

N = 30000;
M = 2;

grand_mean = (mean_chain_1 + mean_chain_2 )/2;

%chain_1 and chain_2 are the chains.

B = N/(M-1) * ((mean_chain_1 - grand_mean).^2 + (mean_chain_2 - grand_mean).^2) ;

v_1 = (std(chain_1)).^2 ;
v_2 = (std(chain_2)).^2 ;

W = (1/M) * (v_1 + v_2) ;



R =  ( (N-1)/N * W + (B/N) ) ./W;


%% plot bars


figure(1)
subplot(2,1,1)
b = bar(R,'FaceColor','flat');

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
line([0.2, 38], [1.1, 1.1], 'Color', [0.1,0.1,0.1],'LineStyle','--',LineWidth=2);
set(gca,'FontSize',24);
xlim([0 38])
ylabel('R_{GR}')
ylim([0 1.2])
yticks([0 0.5 1 1.1]);

set(gca, 'XTick', [5 14 23 30 35]);
set(gca, 'XTickLabel', {'\beta_{ij} ' '\phi_{ij}' '\tau_{ij}' 'r_i'  'D_c'});




%% geweke
store = chainstats(chain_1);
R_geweke = store(:,5);


subplot(2,1,2)
b = bar(R_geweke,'FaceColor','flat');

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
set(gca,'FontSize',24);
xlim([0 38])
ylabel('R_{geweke}')
yticks([0 0.2 0.4 0.6 0.8 0.9 1]);


set(gca, 'XTick', [5 14 23 30 35]);
set(gca, 'XTickLabel', {'\beta_{ij} ' '\phi_{ij}' '\tau_{ij}' 'r_i'  'D_c'});



