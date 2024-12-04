clear all;

addpath(genpath('./../'));

% import data
load('./../community/data/qpcr','data'); % qpcr data
load('./../community/data/parameters_example','pars'); % parameters without nans
pars1 = pars;
load('./../community/data/parameters'); % true parameter set with nans

%% simulate trajectory before inference

model = SEIV_diff_NE_ineffi_infection(5,5,187);
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
model = SEIV_diff_NE_ineffi_infection(5,5,max_NE);
model.name = 'SEIV-noninfectious';

model.host_growth = 0;
model.viral_decay = 0;
model.viral_adsorb = 0;
model.lysis_reset = 0;
model.debris_inhib = 0;




%tvec = 0:0.1:15.75; % for better viz

%% raw data
load('./../community/data/triplicate_data.mat');

tvec = time/60;

hf = figure(1)



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

%clear t1 V1 S1;
max_iter = 10
theta = [5.64e-8 ,1.88e-8 ,   1.24e-7  ,  1e-7 , 3.51e-8,  1.23e-7    5.58e-8, 8.5e-8  ,  4.4e-8,0.19,0.245,0.22,0.28,0.25, 1.9, 126.5  , 63.8  , 36.5, 100,75.1  351.5,87.2  324.1,1.9 ,  1.2500    2.03    2 ,1.8 , 1.45    2.16,  1.37    2.02, 0.2, 0.1, 0.2, 0.1,0.4,0.1, 0.1, 0.1,0.1];
theta_step(1:9) = 1e-9; %phi
theta_step(10:14) = 0.01; %r
theta_step(15:23) = 3; %beta
theta_step(24:32) = 0.1; %tau
theta_step(33:41) = 0.05; %prob

for step = 1:max_iter




pars1.phi = [ 0   theta(1)        0         0         0;
             theta(2)   theta(3)   theta(4)         0         0;
                0         0    theta(5)        0         0;
                0         0         0    theta(6)  theta(7);
                0         0         0   theta(8)    theta(9)];



pars1.r = [theta(10),theta(11),theta(12),theta(13),theta(14)]';



% rewrite pars1
pars1.beta =[0    theta(15)        0         0         0;
  theta(16)   theta(17)   theta(18)       0         0;
         0         0   theta(19)        0         0;
         0         0         0   theta(20)  theta(21);
         0         0         0   theta(22)  theta(23)];



pars1.tau = [0    theta(24)        0         0         0;
  theta(25)   theta(26)   theta(27)       0         0;
         0         0   theta(28)        0         0;
         0         0         0   theta(29)  theta(30);
         0         0         0   theta(31)  theta(32)];



pars1.eta(pars1.tau>0) = 1./pars1.tau(pars1.tau>0);


pars1.prob_infec  = [0    theta(33)        0         0         0;
  theta(34)   theta(35)   theta(36)       0         0;
         0         0   theta(37)        0         0;
         0         0         0   theta(38)  theta(39);
         0         0         0   theta(40)  theta(41)];


pars = pars1;


[t1,S1,V1,~] = simulate_ode_noninfec(model,pars1,tvec,pars1.S0,pars1.V0); % initial parameter set

S1(S1<1e-1) = 0;
V1(V1<1e-1) = 0;
er_host = log(S1) - log([mean(host1');mean(host2');mean(host3');mean(host4');mean(host5')]'*1e3);
er_host_t = sum(abs(er_host(~isinf(er_host))));

er_virus = log(S1) - log([mean(virus1');mean(virus2');mean(virus3');mean(virus4');mean(virus5')]'*1e3);
er_virus_t = sum(abs(er_virus(~isinf(er_virus))));

error_total = er_host_t + er_virus_t


theta_up = theta+theta_step;


end


% for i = 1:5
% subplot(2,5,i)
% plot(t1,S1(:,i),'Color','k','linewidth',2);hold on;
% 
% subplot(2,5,i+5)
% plot(t1,V1(:,i),'Color','k','linewidth',2);hold on;
% end



%% loop over parameter space
% 
% num_iter = 1000;
% 
% for iter = 1:num_iter
% 
% 
% display_statement = sprintf('num iter = %d',iter);
% disp(display_statement);
% 
% prob_random = 0.1 + rand(1,9) *0.5; %uniform b/w 0 and 0.5;
% 
% clear t1 V1 S1;
% 
% 
% 
% 
% 
% [t1,S1,V1,~] = simulate_ode_noninfec(model,pars1,tvec,pars1.S0,pars1.V0); % initial parameter set
% color_of_the_lines=[171,193,157]./255;
% 
% for i = 1:5
% subplot(2,5,i)
% %plot(t1,S1(:,i),'Color',color_of_the_lines,'linewidth',2);hold on;
% patchline(t1,S1(:,i),'linestyle','-','edgecolor',color_of_the_lines,'linewidth',2,'edgealpha',0.1) %S is total bacteria
% hold on;
% 
% 
% subplot(2,5,i+5)
% %plot(t1,V1(:,i),'Color',color_of_the_lines,'linewidth',2);hold on;
% patchline(t1,V1(:,i),'linestyle','-','edgecolor',color_of_the_lines,'linewidth',2,'edgealpha',0.1) %V is total viruses
% hold on;
% 
% end
% 
% end
% 
% %legend('Data','Baseline','Ineffective infection');
% 
% 
% %%
% 
% han=axes(hf,'visible','off'); 
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% set(gca,'FontSize',20);
% set(gca,'fontname','times')  % Set it to times
% xlabel("Time (hours)");
% 
