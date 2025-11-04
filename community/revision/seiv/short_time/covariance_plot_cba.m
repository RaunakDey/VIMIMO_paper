close all;
clear all;
clc;



%%%%  I am going to skip this from the paper now


%% load data.


load("cba_seiv_2.mat",'chain_final');
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

%% figures

figure(1)
for i = 1:5
    for j =1:5

pearson = corrcoef(chain(transient_id:end,i),chain(transient_id:end,j));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if i==j
subplot(5,5,(i-1)*5+j)
smoothHistogram(chain(transient_id:end,i),10,color_input);
axis square
set(gca,'FontSize',14)
elseif i>j
 subplot(5,5,(i-1)*5+j)
 plot(chain(transient_id:end,j),chain(transient_id:end,i),'.','Color',color_input);
 set(gca,'FontSize',14)
 axis square
end

    end
end

%%

figure(2)
for i = 1:5
    for j =1:5

ii = i+5;
jj = j;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(5,5,(i-1)*5+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
 axis square
elseif ii>jj
 subplot(5,5,(i-1)*5+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
 axis square
end

    end
end

%%

figure(3)
for i = 1:5
    for j =1:5

ii = i+5;
jj = j+5;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(5,5,(i-1)*5+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
 axis square
elseif ii>jj
 subplot(5,5,(i-1)*5+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
 axis square
end

    end
end

%%

figure(4)
for i = 1:5
    for j =1:5

ii = i+10;
jj = j+0;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(5,5,(i-1)*5+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
 axis square
elseif ii>jj
 subplot(5,5,(i-1)*5+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
 axis square
end

    end
end


%%

figure(5)
for i = 1:5
    for j =1:5

ii = i+10;
jj = j+5;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(5,5,(i-1)*5+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
 axis square
elseif ii>jj
 subplot(5,5,(i-1)*5+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
 axis square
end

    end
end


%%

figure(6)
for i = 1:5
    for j =1:5

ii = i+10;
jj = j+10;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(5,5,(i-1)*5+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
 axis square
elseif ii>jj
 subplot(5,5,(i-1)*5+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
 axis square
end

    end
end





%% 

figure(7)
for i = 1:3
    for j =1:5

ii = i+15;
jj = j+0;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(5,5,(i-1)*5+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
elseif ii>jj
 subplot(5,5,(i-1)*5+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
end

    end
end

%%

figure(8)
for i = 1:3
    for j =1:5

ii = i+15;
jj = j+5;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(5,5,(i-1)*5+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
elseif ii>jj
 subplot(5,5,(i-1)*5+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
end

    end
end

%%
figure(9)
for i = 1:3
    for j =1:5

ii = i+15;
jj = j+10;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(5,5,(i-1)*5+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
elseif ii>jj
 subplot(5,5,(i-1)*5+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
end

    end
end


%%

figure(10)
for i = 1:3
    for j =1:3

ii = i+15;
jj = j+15;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(5,5,(i-1)*5+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
elseif ii>jj
 subplot(5,5,(i-1)*5+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
end

    end
end

%% joined


axis_labels = ["\beta_{21}","\beta_{12}","\beta_{22}","\beta_{23}","\beta_{33}","log\phi_{21}","log\phi_{12}","log\phi_{22}","log\phi_{23}","log\phi_{33}","\tau_{21}","\tau_{12}","\tau_{22}","\tau_{23}","\tau_{33}","r_1","r_2","r_3"];


figure(11)


for i = 1:18
    for j =1:18

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
subplot(18,18,(i-1)*18+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
elseif ii>jj
 subplot(18,18,(i-1)*18+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
end

if j == 1
    ylabel(axis_labels(i));
end

if i == 18
    xlabel(axis_labels(j));
end

    end
end

