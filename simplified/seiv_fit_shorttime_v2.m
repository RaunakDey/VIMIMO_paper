clear all;
clc;

%%

load("./../community/data/qpcr.mat");
load('./../community/data/triplicate_data.mat');

%% parameters

NH = 3;
NV = 3;


pars.NH = NH;
pars.NV = NV;
pars.M = [ 0     1     0    
     1     1     1   
     0     0     1 ];

pars.NE = [ 0    30     0    
     10    80    71     
     0     0     60];

pars.NE_max = max(pars.NE(:));
NE_max = pars.NE_max;

pars.phi = 1.0e-08 * [0    5.760       0    
   1.122     3.22   0.4528        
         0         0    2.471];

pars.beta = [0    100          0      
            498  400    309.1     
         0         0   106.3];

pars.r = [  0.12
    0.183
    0.244];

            
pars.V0 = [428870 2.8689e+05 528000];
pars.S0 = [2510000 5640000 3.0233e+06];

pars.tau = [0    10        0         
    2  6    6        
         0         0    4.0];

pars.eta = zeros(NH,NV);
pars.eta(pars.tau>0) = 1./pars.tau(pars.tau>0);


%% id


% which variable is which
id.S = 1:NH;
id.E = (1:NH*NV*NE_max)+id.S(end);
id.I = (1:NH*NV)+id.E(end);
id.V = (1:NV)+id.I(end);
id.D = 1+id.V(end); % debris
id.Emat = reshape(id.E,[NH NV NE_max]);
id.Imat = reshape(id.I,[NH NV]);


%% ODE

function dydt = seiv(t,y,pars,id)

            NH = pars.NH;
            NV = pars.NV;
            OH = ones(NH,1);
            OV = ones(NV,1);
            NE_max = pars.NE_max;

            
            S = y(id.S);
            Emat = y(id.Emat);
            Imat = y(id.Imat);
            V = y(id.V);
            D = y(id.D);
            N = S+sum(Emat,3)*OV+Imat*OV;
            etaeff = pars.eta.*(pars.NE+1);  

            dS = pars.r.*S - S.*((pars.M.*pars.phi)*V);
            dEmat = (pars.M.*pars.phi).*(S*V') - etaeff.*Emat(:,:,1);
            dEmat2 = etaeff.*Emat(:,:,1:end-1) - etaeff.*Emat(:,:,2:end);
            for i = 1:NH
                for j = 1:NV
                    if (pars.NE(i,j) ~= pars.NE_max  && pars.NE(i,j) ~=0) 
                        Emat(i,j,pars.NE(i,j)+1:end) = Emat(i,j,pars.NE(i,j));
                    end
                end
            end
            
            dImat = etaeff.*Emat(:,:,end) - etaeff.*Imat;
            dV =  (pars.beta.*etaeff.*Imat)'* OH - V.*((pars.M.*pars.phi)'*N) ;
            dD = sum(etaeff(:).*Imat(:)); % sum across all pairs for net lysis rate
            
            dydt = [dS; dEmat(:); dEmat2(:); dImat(:); dV; dD];
            
        end
    



%% Initial Conditions
S0 = pars.S0(:);  % column vector
E0 = zeros(pars.NH, pars.NV, pars.NE_max);  % all exposed states 0
I0 = zeros(pars.NH, pars.NV);         % all infected states 0
V0 = pars.V0(:);                      % column vector
D0 = 0;
% Flatten y0 into a single column vector
y0 = [S0; E0(:); I0(:); V0; D0];



%% simulate
tvec = 0:0.1:15.75;

% ODE wrapper to include parameters
odefun = @(t, y) seiv(t, y, pars,id);
options = odeset('NonNegative',1,'RelTol',1e-8);
% Run the ODE solver
[t_out, y_out] = ode45(odefun, tvec, y0,options);



%total bacteria
S = y_out(:,id.S);

for i = 1:pars.NH 
    B_out(:,i) = S(:,i) + sum(y_out(:,id.Emat(i,:)),2) + sum(y_out(:,id.Imat(i,:)),2);
end


% virus

V_out = y_out(:,id.V);


%% debug
 model2 =  SEIVD_diff_NE_diff_debris_abs(3,3,80);
 pars.epsilon= ones(1,10);
[t2,B_out,V_out,D_median,I_median,E_median,y] =  simulate_ode(model2,pars,tvec,pars.S0,pars.V0); % mcmc parameter set



%% plot


clear ans;
color_ofthe_fit = [1 0 0]*0.5;
linewidth = 2;

hf4 = figure;
subplot(2,5,1)
errorbar(time/60,mean(1e3*host1'),std(1e3*host1'),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e9]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
ylabel({'Host density';'(cells/ml)'});
title('CBA 4','FontSize',18);
    plot(t_out,B_out(:,1),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,2)
errorbar(time/60,mean(1e3*host2'),std(1e3*host2'),'o','MarkerSize',8,  'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255] );hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e9]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('CBA 18','FontSize',18);
   
    plot(t_out,B_out(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,3)
errorbar(time/60,mean(1e3*host3'),std(1e3*host3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e9]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
title('CBA 38','FontSize',18);

    plot(t_out,B_out(:,3),'-','Color',color_ofthe_fit,'LineWidth',linewidth);








subplot(2,5,6)
errorbar(time/60,mean(1e3*virus1'),std(1e3*virus1'),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e4 1e6 1e8 1e10]);
ylabel({'Phage density';'(virions/ml)'});
title('\phi18:2','FontSize',18);
    
    plot(t_out,V_out(:,1),'-','Color',color_ofthe_fit,'LineWidth',linewidth);

subplot(2,5,7)
errorbar(time/60,mean(1e3*virus2'),std(1e3*virus2'),'o','MarkerSize',8,  'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255] );hold on;
set(gca, 'YScale', 'log');
set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
    yticks([1e4 1e6 1e8 1e10]);
    title('\phi18:3','FontSize',18);
    
    plot(t_out,V_out(:,2),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,8)
errorbar(time/60,mean(1e3*virus3'),std(1e3*virus3'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('\phi38:1','FontSize',18);
    
    plot(t_out,V_out(:,3),'-','Color',color_ofthe_fit,'LineWidth',linewidth);


han=axes(hf4,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
xlabel("Time (hours)");

