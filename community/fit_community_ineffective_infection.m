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
    
%%

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
    1.2500    2.03    2         0         0;
         0         0    1.8         0         0;
         0         0         0    1.45    2.16 ;
         0         0         0    1.37    2.02];

pars1.eta(pars1.tau>0) = 1./pars1.tau(pars1.tau>0);


pars = pars1;
% no inhibition
pars1.prob_infec = 1*[0 1 0 0 0;
                        1 1 1 0 0;
                        0 0 1 0 0;
                        0 0 0 1 1;
                        0 0 0 1 1];

[t1,S1,V1,~] = simulate_ode_noninfec(model,pars1,tvec,pars1.S0,pars1.V0); % initial parameter set


for i = 1:5
subplot(2,5,i)
plot(t1,S1(:,i),'Color','k','linewidth',2);hold on;

subplot(2,5,i+5)
plot(t1,V1(:,i),'Color','k','linewidth',2);hold on;
end

%% figures

%clear t1 V1 S1;
max_iter = 1;

%theta = [5.64e-8 ,1.88e-8 ,   1.24e-7  ,  1e-7 , 3.51e-8,  1.23e-7    5.58e-8, 8.5e-8  ,  4.4e-8,0.19,0.245,0.22,0.28,0.25, 1.9, 126.5  , 63.8  , 36.5, 100,75.1  351.5,87.2  324.1,1.9 ,  1.2500    2.03    2 ,1.8 , 1.45    2.16,  1.37    2.02, 0.2, 0.1, 0.2, 0.1,0.4,0.1, 0.1, 0.1,0.1];

%theta = [5.50000000000000e-08	1.74000000000000e-08	1.22600000000000e-07	9.85999999999997e-08	3.37000000000000e-08	1.21600000000000e-07	5.44000000000000e-08	8.35999999999997e-08	4.26000000000000e-08	0.189999998600000	0.244999998600000	0.219999998600000	0.279999998600000	0.249999998600000	1.89999999860000	126.499999998600	63.7999999986001	36.4999999986001	99.9999999985998	75.0999999985998	351.499999998601	87.1999999985998	324.099999998601	1.89999999860000	1.24999999860000	2.02999999860000	1.99999999860000	1.79999999860000	1.44999999860000	2.15999999860000	1.36999999860000	2.01999999860000	0.199999998600000	0.0999999985999999	0.199999998600000	0.0999999985999999	0.399999998600000	0.0999999985999999	0.0999999985999999	0.0999999985999999	0.0999999985999999]
theta = [5.49200000000000e-08	1.73200000000000e-08	1.22520000000000e-07	9.85199999999997e-08	3.36200000000000e-08	1.21520000000000e-07	5.43200000000000e-08	8.35199999999997e-08	4.25200000000000e-08	0.189199998600000	0.244199998600000	0.219199998600000	0.279199998600000	0.249199998600000	1.65999999860000	126.259999998600	63.5599999986001	36.2599999986001	99.7599999985998	74.8599999985998	351.259999998601	86.9599999985998	323.859999998601	1.89199999860000	1.24199999860000	2.02199999860000	1.99199999860000	1.79199999860000	1.44199999860000	2.15199999860000	1.36199999860000	2.01199999860000	0.195999998600000	0.0959999985999999	0.195999998600000	0.0959999985999999	0.395999998600000	0.0959999985999999	0.0959999985999999	0.0959999985999999	0.0959999985999999];

theta_step(1:9) = 1e-9; %phi
theta_step(10:14) = 0.025; %r
theta_step(15:23) = 3; %beta
theta_step(24:32) = 0.1; %tau
theta_step(33:41) = 0.025; %prob



theta_upper(1:9) = 1e-7*ones(1,9);
theta_lower(1:9) = 1e-10*ones(1,9);
theta_upper(10:14) = 0.7*ones(1,5); %r
theta_lower(10:14) = 0.1*ones(1,5); %r
theta_upper(15:23) = ones(1,9)*600;
theta_lower(15:23) = ones(1,9)*10;
theta_upper(24:32) = 10*ones(1,9); %tau
theta_lower(24:32) = 1*ones(1,9);
theta_upper(33:41) = ones(1,9); 
theta_lower(33:41) = 0.1*ones(1,9);





% theta_step(1:9) = 1e-9; %phi
% theta_step(10:14) = 0.05; %r
% theta_step(15:23) = 3; %beta
% theta_step(24:32) = 0.2; %tau
% theta_step(33:41) = 0.05; %prob




learning_rate = 0.1;  % Initial learning rate
tolerance = 1e-6;     % Tolerance for stopping criterion
max_iter = 15;       % Maximum number of iterations


for j = 1:10

theta = theta_lower + rand(1,41).*(theta_upper -theta_lower) ;

for step = 1:max_iter
    
    step
    % Calculate the current error
    er_t = error_ineffi_seiv(theta, pars1, model, tvec, host1, host2, host3, host4, host5, virus1, virus2, virus3, virus4, virus5)
    
    % for j=1:41
    %     % Calculate theta updates for finite differences (numerical gradient)
    %     theta_up = theta;
    %     theta_down = theta;
    % 
    %     theta_up(j) = theta(j) +learning_rate * theta_step(j);
    %     theta_down(j) = theta(j) - learning_rate * theta_step(j);
    % 
    %     er_t_up = error_ineffi_seiv(theta_up, pars1, model, tvec, host1, host2, host3, host4, host5, virus1, virus2, virus3, virus4, virus5);
    %     er_t_down = error_ineffi_seiv(theta_down, pars1, model, tvec, host1, host2, host3, host4, host5, virus1, virus2, virus3, virus4, virus5);
    %     del_C_del_theta(j) = (er_t_up - er_t_down) ./ (theta_up(j) - theta_down(j));
    %     j
    % end

    
 
     theta_up = theta +learning_rate * theta_step;
     theta_down = theta - learning_rate * theta_step;

     er_t_up = error_ineffi_seiv(theta_up, pars1, model, tvec, host1, host2, host3, host4, host5, virus1, virus2, virus3, virus4, virus5);
     er_t_down = error_ineffi_seiv(theta_down, pars1, model, tvec, host1, host2, host3, host4, host5, virus1, virus2, virus3, virus4, virus5);
  
     grad = (er_t_up - er_t_down) ./ (theta_up - theta_down);
     grad_normalized = grad./sqrt(sum(grad.^2));

     signature = sign(er_t_up-er_t_down);

    % Update theta using gradient descent rule
    %theta = theta - learning_rate * grad_normalized.*theta_step;
    theta = theta - learning_rate * signature.*theta_step;
    %theta = theta - learning_rate*del_C_del_theta_normalized;
    
    
    % Check for convergence (if gradient is small enough)
    % if norm(grad) < tolerance
    %     fprintf('Converged at iteration %d\n', step);
    %     break;
    % end
    
end

theta_final(:,j) = theta;
error(j) =  er_t;

end



%%
% Output final theta values
[value,index] = min(error);
theta = theta_final(:,index);


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




clear t1 V1 S1;



[t1,S1,V1,~] = simulate_ode_noninfec(model,pars1,tvec,pars1.S0,pars1.V0); % initial parameter set
color_of_the_lines=[171,193,157]./255;

for i = 1:5
subplot(2,5,i)
%plot(t1,S1(:,i),'Color',color_of_the_lines,'linewidth',2);hold on;
patchline(t1,S1(:,i),'linestyle','-','edgecolor',color_of_the_lines,'linewidth',2,'edgealpha',0.9) %S is total bacteria
hold on;


subplot(2,5,i+5)
%plot(t1,V1(:,i),'Color',color_of_the_lines,'linewidth',2);hold on;
patchline(t1,V1(:,i),'linestyle','-','edgecolor',color_of_the_lines,'linewidth',2,'edgealpha',0.9) %V is total viruses
hold on;


end



han=axes(hf,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
xlabel("Time (hours)");

