clear all;
clc;

skip  = 1 ; % from chain thinning.
% if not in use set to 1.

load("v27.mat",'chain');
chain_1 = chain(1:skip:end,:);



load("v25.mat",'chain');
chain_2 = chain(1:skip:end,:);

fs = 22;

% 
% load("v25.mat",'chain');
% chain_1 = chain(30001:40000,:);
% 
% chain_2 = chain(40001:end,:);
%% for beta

fig1 = figure;

transient_id =2000;% 20000;
 

color1 = [76,132,147]./255;
color2 = [217,76,33]./255;

for i = 1:9

subplot(5,2,i)
plot(chain_1(transient_id:end,i),'Color',color1 );

hold on;
set(gca,'FontSize',fs)
end


for i = 1:9

subplot(5,2,i)
plot(chain_2(transient_id:end,i),'Color',color2);
%ylim([  min(min(chain_1(transient_id:end,5+i)),min(chain_2(transient_id:end,5+i)))-0.1,  max(max(chain_1(transient_id:end,5+i)),max(chain_2(transient_id:end,5+i)))+0.1 ]);

make_title(i);
set(gca,'FontSize',fs)
xlim([0 length(chain_2(transient_id:end,i))])
end

han=axes(fig1,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'\beta (virions/cell)');
xlabel(han,'steps');
%title(han,'yourTitle');
set(gca,'FontSize',fs)
%% for r and Dc



fig2 = figure;

 
for i = 1:5

subplot(5,2,2*i-1)
plot(chain_1(transient_id:end,i+27), 'Color',color1 );

hold on;
set(gca,'FontSize',fs)
end


for i = 1:5

subplot(5,2,2*i-1)
plot(chain_2(transient_id:end,i+27),'Color',color2);
%ylim([  min(min(chain_1(transient_id:end,i)),min(chain_2(transient_id:end,i))),  max(max(chain_1(transient_id:end,i)),max(chain_2(transient_id:end,i))) ]);

make_title_host(i);
xlim([0 length(chain_2(transient_id:end,i))])

if i ==3
    ylabel('r (cells/hr)')
end
set(gca,'FontSize',fs)


end

for i = 1:5

subplot(5,2,2*i)
plot(chain_1(transient_id:end,i+32), 'Color',color1 );

hold on;
set(gca,'FontSize',fs)
end


for i = 1:5

subplot(5,2,2*i)
plot(chain_2(transient_id:end,i+32),'Color',color2);
%ylim([  min(min(chain_1(transient_id:end,i)),min(chain_2(transient_id:end,i))),  max(max(chain_1(transient_id:end,i)),max(chain_2(transient_id:end,i))) ]);

make_title_host(i);
xlim([0 length(chain_2(transient_id:end,i))])


if i ==3
    ylabel('log_{10}D_c (dead cells/ml)')
end

set(gca,'FontSize',fs)
end





han=axes(fig2,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'r (cells/hr)');
xlabel(han,'steps');
%title(han,'yourTitle');
set(gca,'FontSize',fs)

%% for phi

fig3 = figure;


for i = 1:9
subplot(5,2,i)
plot(chain_1(transient_id:end,9+i), 'Color',color1 );
hold on;
set(gca,'FontSize',fs)
end



for i = 1:9
subplot(5,2,i)
plot(chain_2(transient_id:end,9+i),'Color',color2);
%ylim([  min(min(chain_1(transient_id:end,14+i)),min(chain_2(transient_id:end,14+i)))-0.1,  max(max(chain_1(transient_id:end,14+i)),max(chain_2(transient_id:end,14+i)))+0.1 ]);
make_title(i);
xlim([0 length(chain_2(transient_id:end,i))])
set(gca,'FontSize',fs)
end

han=axes(fig3,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'log_{10} \phi (ml/hr) ');
xlabel(han,'steps');
%title(han,'yourTitle');
set(gca,'FontSize',fs)


%% for tau


fig4 = figure;

for i = 1:9

subplot(5,2,i)
plot(chain_1(transient_id:end,18+i),'Color',color1);

hold on;
set(gca,'FontSize',fs)
end



for i = 1:9

subplot(5,2,i)
plot(chain_2(transient_id:end,18+i),'Color',color2);
%ylim([  min(min(chain_1(transient_id:end,24+i)),min(chain_2(transient_id:end,24+i)) -0.1) ,  max(max(chain_1(transient_id:end,24+i)),max(chain_2(transient_id:end,24+i)) +0.1)]);
xlim([0 length(chain_2(transient_id:end,i))])
set(gca,'FontSize',fs)
make_title(i);
end


han=axes(fig4,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'\tau (hr) ');
xlabel(han,'steps');
set(gca,'FontSize',fs)
%title(han,'yourTitle');
