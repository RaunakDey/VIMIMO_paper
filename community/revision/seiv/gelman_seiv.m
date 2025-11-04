clear all;
clc;

skip  = 5 ; % from chain thinning.
% if not in use set to 1.

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


chain_1 = chain_joined(20000:skip:35000,:);
chain_2 = chain_joined(35000:skip:50000,:);

mean_chain_1 = mean(chain_1);
mean_chain_2 = mean(chain_2);



%%

N = length(chain_1);
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





set(gca,'XTickLabel','');
hold on;
line([0.2, 38], [1.1, 1.1], 'Color', [0.1,0.1,0.1],'LineStyle','--',LineWidth=2);
set(gca,'FontSize',24);
xlim([0 34])
ylabel('R_{GR}')
ylim([0 1.2])
yticks([0 0.5 1 1.1]);

set(gca, 'XTick', [5 14 23 30]);
set(gca, 'XTickLabel', {'\beta_{ij} ' '\phi_{ij}' '\tau_{ij}' 'r_i'});




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



set(gca,'XTickLabel','');
hold on;
line([0.2, 38], [0.9, 0.9],'Color', [0.1,0.1,0.1],'LineStyle','--',LineWidth=2  );
set(gca,'FontSize',24);
xlim([0 34])
ylabel('R_{geweke}')
yticks([0 0.2 0.4 0.6 0.8 0.9 1]);


set(gca, 'XTick', [5 14 23 30]);
set(gca, 'XTickLabel', {'\beta_{ij} ' '\phi_{ij}' '\tau_{ij}' 'r_i' });



