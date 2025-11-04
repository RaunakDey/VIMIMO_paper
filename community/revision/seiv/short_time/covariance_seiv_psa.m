close all;
clear all;
clc;



%%%%  I am going to skip this from the paper now


%% load data.


load("psa_seiv_2.mat",'chain_final');
chain_cba = chain_final; 


skip = 5;
% if not used please set skip to 1


transient_id = 25000/skip;


chain = chain_cba(1:skip:end,:);




%% color plots

x = -1:0.01:1;
y= ones(length(x),1);
for i = 1:length(x)
pearson_corr = x(i);

if x(i)>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
elseif x(i)<0
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

scatter(x(i),y(i),MarkerFaceColor=color_input,MarkerEdgeColor=color_input);
hold on;

end
set(gca,'box','off')

yticks([]);
set(gca,'FontSize',14)

%%

axis_labels = ["\beta_{44}","\beta_{54}","\beta_{45}","\beta_{55}","log\phi_{44}","log\phi_{54}","log\phi_{45}","log\phi_{55}","\tau_{44}","\tau_{54}","\tau_{45}","\tau_{55}","r_4","r_5"];



figure(11)


for i = 1:14
    for j =1:14

ii = i+0;
jj = j+0;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(14,14,(i-1)*14+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
elseif ii>jj
 subplot(14,14,(i-1)*14+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
end

if j == 1
    ylabel(axis_labels(i));
end

if i == 14
    xlabel(axis_labels(j));
end

    end
end

