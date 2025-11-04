close all;
clear all;
clc;

%% load data.

load('./v25.mat');
close all;

skip = 10;
% if not used plkease set skip to 1


transient_id = 25000/skip;
chain = chain(1:skip:end,:);

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
for i = 1:9
    for j =1:9

pearson = corrcoef(chain(transient_id:end,i),chain(transient_id:end,j));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if i==j
subplot(9,9,(i-1)*9+j)
smoothHistogram(chain(transient_id:end,i),10,color_input);
axis square
set(gca,'FontSize',14)
elseif i>j
 subplot(9,9,(i-1)*9+j)
 plot(chain(transient_id:end,j),chain(transient_id:end,i),'.','Color',color_input);
 set(gca,'FontSize',14)
 axis square
end

    end
end


figure(2)
for i = 1:9
    for j =1:9

ii = i+9;
jj = j;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(9,9,(i-1)*9+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
 axis square
elseif ii>jj
 subplot(9,9,(i-1)*9+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
 axis square
end

    end
end


figure(3)
for i = 1:9
    for j =1:9

ii = i+9;
jj = j+9;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(9,9,(i-1)*9+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
 axis square
elseif ii>jj
 subplot(9,9,(i-1)*9+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
 axis square
end

    end
end


figure(4)
for i = 1:9
    for j =1:9

ii = i+18;
jj = j+0;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(9,9,(i-1)*9+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
 axis square
elseif ii>jj
 subplot(9,9,(i-1)*9+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
 axis square
end

    end
end




figure(5)
for i = 1:9
    for j =1:9

ii = i+18;
jj = j+9;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(9,9,(i-1)*9+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
 axis square
elseif ii>jj
 subplot(9,9,(i-1)*9+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
 axis square
end

    end
end




figure(6)
for i = 1:9
    for j =1:9

ii = i+18;
jj = j+18;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(9,9,(i-1)*9+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
 axis square
elseif ii>jj
 subplot(9,9,(i-1)*9+j)
 plot(chain(transient_id:end,jj),chain(transient_id:end,ii),'.','Color',color_input);set(gca,'FontSize',14)
 axis square
end

    end
end







%% 

figure(7)
for i = 1:9
    for j =1:5

ii = i+27;
jj = j+0;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(9,9,(i-1)*9+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
elseif ii>jj
 subplot(9,9,(i-1)*9+j)
 plot(chain(transient_id:end,ii),chain(transient_id:end,jj),'.','Color',color_input);set(gca,'FontSize',14)
end

    end
end


figure(8)
for i = 1:9
    for j =1:5

ii = i+27;
jj = j+9;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(9,9,(i-1)*9+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
elseif ii>jj
 subplot(9,9,(i-1)*9+j)
 plot(chain(transient_id:end,ii),chain(transient_id:end,jj),'.','Color',color_input);set(gca,'FontSize',14)
end

    end
end


figure(9)
for i = 1:9
    for j =1:5

ii = i+27;
jj = j+18;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(9,9,(i-1)*9+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
elseif ii>jj
 subplot(9,9,(i-1)*9+j)
 plot(chain(transient_id:end,ii),chain(transient_id:end,jj),'.','Color',color_input);set(gca,'FontSize',14)
end

    end
end


figure(10)
for i = 1:5
    for j =1:5

ii = i+27;
jj = j+27;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(9,9,(i-1)*9+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
elseif ii>jj
 subplot(9,9,(i-1)*9+j)
 plot(chain(transient_id:end,ii),chain(transient_id:end,jj),'.','Color',color_input);set(gca,'FontSize',14)
end

    end
end



figure(15)
for i = 1:5
    for j =1:5

ii = i+32;
jj = j+32;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(9,9,(i-1)*9+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
elseif ii>jj
 subplot(9,9,(i-1)*9+j)
 plot(chain(transient_id:end,ii),chain(transient_id:end,jj),'.','Color',color_input);set(gca,'FontSize',14)
end

    end
end



%%

figure(11)
for i = 1:5
    for j =1:9

ii = i+32;
jj = j+0;
pearson = corrcoef(chain(transient_id:end,ii),chain(transient_id:end,jj));
pearson_corr = pearson(1,2);

if pearson_corr>0
color_input = pearson_corr*[1 0 0]   + (1-pearson_corr)*[1 1 0];
else
color_input = (-pearson_corr)*[0 1 0] + (1+pearson_corr)*[1 1 0] ;
end

if ii==jj
subplot(9,9,(i-1)*9+j)
smoothHistogram(chain(transient_id:end,ii),10,color_input);set(gca,'FontSize',14)
elseif ii>jj
 subplot(9,9,(i-1)*9+j)
 plot(chain(transient_id:end,ii),chain(transient_id:end,jj),'.','Color',color_input);set(gca,'FontSize',14)
end

    end



end