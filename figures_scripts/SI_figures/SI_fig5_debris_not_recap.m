%%%%% later.

clear all;

addpath(genpath('./../'));

% import data
load('./../../community/data/qpcr','data'); % qpcr data
load('./../../community/data/parameters_example','pars'); % parameters without nans
pars1 = pars;
load('./../../community/data/parameters'); % true parameter set with nans

%% simulate trajectory before inference



flags.phi_entire_matrix = 0;
flags.ssfun_normalized = 0;
flags.tau_mult = 1;
flags.mcmc_algorithm = 1; % default is 1 ('dram')
flags.inference_script = 1;
flags.confidence_interval = 0;
flags.tau_new = 0;

pars.NE = 10*(pars.M == 1);
pars.NE(2,1) = 135;
pars.NE(1,2) = 10;
pars.NE(2,2) = 172;
pars.NE(2,3) = 100;
pars.NE(3,3) = 77;
pars.NE(4,4) = 138;
pars.NE(5,4) = 108;
pars.NE(4,5) = 167;
pars.NE(5,5) = 187;
pars1.NE = pars.NE;
max_NE = round(max(max(pars.NE)));
model =  SEIVD_diff_NE_diff_debris_abs(5,5,max_NE);
model.host_growth = 0;
model.viral_decay = 0;
model.viral_adsorb = 0;
model.lysis_reset = 0;
model.debris_inhib = 1;
model.debris_inhib2 = 1;
model.debris_inhib3 = 1;
model.debris_inhib4 = 1;
model.debris_inhib5 = 1;


load('./../../community/results/v25.mat','pars2');

pars1.Dc = pars2.Dc;
pars1.Dc2 = pars2.Dc2;
pars1.Dc3 = pars2.Dc3;
pars1.Dc4 = pars2.Dc4;
pars1.Dc5 = pars2.Dc5;


tvec = 0:0.1:15.75; % for better viz

%% raw data
load('./../../community/data/triplicate_data.mat');

hf = figure



  subplot(2,5,1)
errorbar(time/60,mean(1e3*host1'),std(1e3*host1'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255], Color=[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 4','FontSize',18);
    
  

subplot(2,5,2)
errorbar(time/60,mean(1e3*host2'),std(1e3*host2'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255], Color=[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 18','FontSize',18);
    


subplot(2,5,3)
errorbar(time/60,mean(1e3*host3'),std(1e3*host3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255], Color=[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 38','FontSize',18);
    

  

subplot(2,5,6)
errorbar(time/60,mean(1e3*virus1'),std(1e3*virus1'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255], Color=[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('\phi18:1','FontSize',18);
   


subplot(2,5,7)
errorbar(time/60,mean(1e3*virus2'),std(1e3*virus2'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255], Color=[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('\phi18:3','FontSize',18);
    

subplot(2,5,8)
errorbar(time/60,mean(1e3*virus3'),std(1e3*virus3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255], Color=[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('\phi38:1','FontSize',18);



subplot(2,5,4)
errorbar(time/60,mean(1e3*host4'),std(1e3*host4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],Color = [70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA H100','FontSize',18);
    
  

subplot(2,5,5)
errorbar(time/60,mean(1e3*host5'),std(1e3*host5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],Color = [70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e8]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA 13-15','FontSize',18);
    
  

subplot(2,5,9)
errorbar(time/60,mean(1e3*virus4'),std(1e3*virus4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],Color = [70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('PSA HP1','FontSize',18);
   


subplot(2,5,10)
errorbar(time/60,mean(1e3*virus5'),std(1e3*virus5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],Color = [70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('PSA HS6','FontSize',18);
    

%% figures


%%%% these are the correct values.
pars1.phi = [ 0    5.64e-8         0         0         0;
             1.88e-8    1.24e-7    1e-7         0         0;
                0         0    3.51e-8         0         0;
                0         0         0    1.23e-7    5.58e-8;
                0         0         0    8.5e-8    4.4e-8];



pars1.r = [0.19,0.245,0.22,0.28,0.25]';


% rewrite pars1
pars1.beta =[0    1.9         0         0         0;
   126.5   63.8    100         0         0;
         0         0   36.5         0         0;
         0         0         0   75.1  351.5;
         0         0         0   87.2  324.1];



pars1.tau = [0    1.9         0         0         0;
    1.2500    2.03    1        0         0;
         0         0    1.8         0         0;
         0         0         0    1.45    2.16 ;
         0         0         0    1.37    2.02];




pars1.eta(pars1.tau>0) = 1./pars1.tau(pars1.tau>0);


pars_sample = pars1;


for count = 1:500
count
pars_sample.Dc = 10^(log(pars1.Dc)/log(10) + (rand - 0.5)); % in log space \sigma = 0.5;
pars_sample.Dc2 = 10^(log(pars1.Dc2)/log(10) + (rand - 0.5));
pars_sample.Dc3 = 10^(log(pars1.Dc3)/log(10) + (rand - 0.5));
pars_sample.Dc4 = 10^(log(pars1.Dc4)/log(10) + (rand - 0.5));
pars_sample.Dc5 = 10^(log(pars1.Dc5)/log(10) + (rand - 0.5));

[t1,S1,V1,~] = simulate_ode(model,pars_sample,tvec,pars_sample.S0,pars_sample.V0); % initial parameter set

color_of_the_lines=[250 218 221]/255;

for i = 1:5
subplot(2,5,i)
patchline(t1,S1(:,i),'edgecolor',color_of_the_lines,'linewidth',3,'edgealpha',0.5);hold on;

subplot(2,5,i+5)
patchline(t1,V1(:,i),'edgecolor',color_of_the_lines,'linewidth',3,'edgealpha',0.5);hold on;
end

end




han=axes(hf,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
xlabel("Time (hours)");

