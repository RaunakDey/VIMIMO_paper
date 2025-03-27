clc;    
clear all;

load('./../community/results/v25.mat');
load('./../community/data/triplicate_data.mat');
%color_of_the_fit = [250 218 221]/255;
color_of_the_fit = [171,193,157]./255;
blue_color_old = [70/255,130/255,180/255];

%%
transient_id = 30000;

close all;

red1 = [255,186,186]/255;
red2 = [255,123,123]/255;
red3 = [255,82,82]/255;
red4 = [255,0,0]/255;
red5 = [167,0,0]/255;


blue1 = [179,205,224]/255;
blue2 = [100,151,177]/255;
blue3 = [0,91,150]/255;
blue4 = [3,10,208]/255;
blue5 = [0,0,75]/255;

transparency = 0.2;
color = [70/255,130/255,180/255];

figure(1)

%%
subplot(2,5,1)
errorbar(time/60,mean(1e3*host1'),std(1e3*host1'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue1, Color=blue1);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 4','FontSize',18);
    
  

subplot(2,5,2)
errorbar(time/60,mean(1e3*host2'),std(1e3*host2'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue2, Color=blue2);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 18','FontSize',18);
    


subplot(2,5,3)
errorbar(time/60,mean(1e3*host3'),std(1e3*host3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue3, Color=blue3);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 38','FontSize',18);
    

  




subplot(2,5,4)
errorbar(time/60,mean(1e3*host4'),std(1e3*host4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue4,Color = blue4);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA H100','FontSize',18);
    
  

subplot(2,5,5)
errorbar(time/60,mean(1e3*host5'),std(1e3*host5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',blue5,Color = blue5);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA 13-15','FontSize',18);
    
  
    subplot(2,5,6)
errorbar(time/60,mean(1e3*virus1'),std(1e3*virus1'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red1, Color=red1);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('\phi18:1','FontSize',18);
   


subplot(2,5,7)
errorbar(time/60,mean(1e3*virus2'),std(1e3*virus2'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red2, Color=red2);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('\phi18:3','FontSize',18);
    

subplot(2,5,8)
errorbar(time/60,mean(1e3*virus3'),std(1e3*virus3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red3, Color=red3);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('\phi38:1','FontSize',18);


subplot(2,5,9)
errorbar(time/60,mean(1e3*virus4'),std(1e3*virus4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red4,Color = red4);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('PSA HP1','FontSize',18);
   


subplot(2,5,10)
errorbar(time/60,mean(1e3*virus5'),std(1e3*virus5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',red5,Color = red5);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('PSA HS6','FontSize',18);

  %%

  % that 95% confidence interval limit.
  % you can run for all traces and find the quantile at which that
  % "Half-Gaussian" is cutoff (half-gaussian is an approximation).

  cutoff = 200;
   

count=1;

%%

for i = 30000:100:50000
    i
    
    error_store(count) = loglikefun(chain(i,:),data,pars2,mcmcpars,model,0);
    count=count+1;

    [ss, err_vec,~,~,~] = loglikefun(chain(i,:),data,pars2,mcmcpars,model,0);
if loglikefun(chain(i,:),data,pars2,mcmcpars,model,0) < cutoff

 pars_samples = update_pars(pars2,chain(i,:),mcmcpars);
 pars_samples.beta2 = pars_samples.beta;

[t3,S3,V3,~] = simulate_ode(model,pars_samples,tvec,pars2.S0,pars2.V0);
    figure(1)
    subplot(2,5,6)
    patchline(t3,V3(:,1),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;
  
    subplot(2,5,7)
    patchline(t3,V3(:,2),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;

    subplot(2,5,8)
    patchline(t3,V3(:,3),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;

    subplot(2,5,9)
    patchline(t3,V3(:,4),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;

    subplot(2,5,10)
    patchline(t3,V3(:,5),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;

    
    subplot(2,5,1)
    patchline(t3,S3(:,1),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;

    subplot(2,5,2)
    patchline(t3,S3(:,2),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;
    
    subplot(2,5,3)
    patchline(t3,S3(:,3),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;
    
    subplot(2,5,4)
    patchline(t3,S3(:,4),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;
    
    subplot(2,5,5)
    patchline(t3,S3(:,5),'edgecolor',color_of_the_fit,'linewidth',5,'edgealpha',transparency);hold on;

end

end



%% touch up

figure(1)

subplot(2,5,1)
hold on;
set(gca,'FontSize',18);
title('CBA 4');
ylabel({'Host density','(cells/ml)'})
xticks(0:2:16);
ylim([1e4 1e9]);
yticks([1e4 1e5, 1e6, 1e7, 1e8 1e9]);
xlim([0 16]);
axis square;

subplot(2,5,2)
hold on;
set(gca,'FontSize',18);
title('CBA 18');
xticks(0:2:16);
ylim([1e4 1e9]);
yticks([1e4 1e5, 1e6, 1e7, 1e8 1e9]);
xlim([0 16]);
axis square;

subplot(2,5,3)
hold on;
set(gca,'FontSize',18);
title('CBA 38');
xticks(0:2:16);
ylim([1e4 1e9]);
yticks([1e4 1e5, 1e6, 1e7, 1e8 1e9]);
xlim([0 16]);
axis square;


subplot(2,5,4)
hold on;
set(gca,'FontSize',18);
title('PSA H100');
xticks(0:2:16);
ylim([1e4 1e9]);
yticks([1e4 1e5, 1e6, 1e7, 1e8 1e9]);
xlim([0 16]);
axis square;

subplot(2,5,5)
hold on;
set(gca,'FontSize',18);
title('PSA 13-15');
xticks(0:2:16);
ylim([1e4 1e9]);
yticks([1e4 1e5, 1e6, 1e7, 1e8 1e9]);
xlim([0 16]);
axis square;


subplot(2,5,6)
hold on;
set(gca,'FontSize',18);
title('\phi18:2');
ylabel({'Phage density','(virions/ml)'})
xticks(0:2:16);
ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]);
xlim([0 16]);
axis square;



subplot(2,5,7)
hold on;
set(gca,'FontSize',18);
title('\phi18:3');
xticks(0:2:16);
ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]);
xlim([0 16]);
axis square;


subplot(2,5,8)
hold on;
set(gca,'FontSize',18);
title('\phi38:1');
xticks(0:2:16);
ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]);
xlabel('Time (hrs)');
xlim([0 16]);
axis square;


subplot(2,5,9)
hold on;
set(gca,'FontSize',18);
title('PSA HP1');
xticks(0:2:16);
ylim([1e4 1e11]);
xlim([0 16]);
yticks([1e4, 1e6, 1e8, 1e10]);
axis square;


subplot(2,5,10)
hold on;
set(gca,'FontSize',18);
title('PSA HS6');
xticks(0:2:16);
ylim([1e4 1e11]);
yticks([1e4, 1e6, 1e8, 1e10]);
xlim([0 16]);
axis square;
