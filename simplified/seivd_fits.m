clear all;
clc;

%%

load("./../community/data/qpcr.mat");
load('./../community/data/triplicate_data.mat');

%% parameters

NH = 5;
NV = 5;
NE_max = 200;

pars.NH = NH;
pars.NV = NV;
pars.M = [ 0     1     0     0     0
     1     1     1     0     0
     0     0     1     0     0
     0     0     0     1     1
     0     0     0     1     1];

pars.NE = NE_max*pars.M;
pars.NE_max = NE_max;

pars.phi = 1.0e-07 * [0    0.5892         0         0         0
    0.1532    0.7788    0.2407         0         0
         0         0    0.7904         0         0
         0         0         0    0.6114    0.1191
         0         0         0    0.6012    0.2264];

pars.beta = [0    2.8875         0         0         0
  194.9000  204.5800  100.4600         0         0
         0         0   19.9220         0         0
         0         0         0  525.3800   60.6530
         0         0         0  488.0600   51.3260];

pars.r = [ 0.1769
    0.2207
    0.2939
    0.6658
    0.5281];


pars.Dc1 = 5041500;
pars.Dc2 =  5962700;
pars.Dc3 = 11927000;
pars.Dc4 =  1832400;
pars.Dc5 = 1520000;
            
pars.V0 = [428870 2.8689e+05 528000 1.1033e+05 11510000];
pars.S0 = [2510000 5640000 3.0233e+06 6.2033e+06 7.7533e+06];

pars.tau = [0    2.9978         0         0         0
    1.7394    2.7460    2.3189         0         0
         0         0    1.9880         0         0
         0         0         0    1.8220    4.7139
         0         0         0    2.3157    1.9868];

pars.eta = [0    0.3333         0         0         0
    0.5882    0.3704    0.4348         0         0
         0         0    0.5000         0         0
         0         0         0    0.5556    0.2128
         0         0         0    0.4348    0.5000];

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

function dydt = seivd(t,y,pars,id)

            NH = pars.NH;
            NV = pars.NV;
            OH = ones(NH,1);
            OV = ones(NV,1);
            NE_max = 200;

            
            S = y(id.S);
            Emat = y(id.Emat);
            Imat = y(id.Imat);
            V = y(id.V);
            D = y(id.D);
            N = S+sum(Emat,3)*OV+Imat*OV;
            etaeff = pars.eta.*(pars.NE_max+1);  


            Sdeb(1) =  y(1)*1./(1+(D./pars.Dc1).^2); % debris
            Sdeb(2) =  y(2)*1./(1+(D./pars.Dc2).^2);
            Sdeb(3) =  y(3)*1./(1+(D./pars.Dc3).^2);
            Sdeb(4) =  y(4)*1./(1+(D./pars.Dc4).^2);
            Sdeb(5) =  y(5)*1./(1+(D./pars.Dc5).^2);
            

            Ndeb(1) =  N(1)*1./(1+(D./pars.Dc1).^2); % debris
            Ndeb(2) =  N(2)*1./(1+(D./pars.Dc2).^2);
            Ndeb(3) =  N(3)*1./(1+(D./pars.Dc3).^2);
            Ndeb(4) =  N(4)*1./(1+(D./pars.Dc4).^2);
            Ndeb(5) =  N(5)*1./(1+(D./pars.Dc5).^2);
            
                       
            
            Sdeb = Sdeb';
            
            dS = pars.r.*S - Sdeb.*((pars.M.*pars.phi)*V);
            dEmat = (pars.M.*pars.phi).*(Sdeb*V') - etaeff.*Emat(:,:,1);
            dEmat2 = etaeff.*Emat(:,:,1:end-1) - etaeff.*Emat(:,:,2:end);
            for i = 1:NH
                for j = 1:NV
                    if (pars.NE(i,j) ~= pars.NE_max  && pars.NE(i,j) ~=0) 
                        Emat(i,j,pars.NE(i,j)+1:end) = Emat(i,j,pars.NE(i,j));
                    end
                end
            end
            
            dImat = etaeff.*Emat(:,:,end) - etaeff.*Imat;
            dV =  (pars.beta.*etaeff.*Imat)'* OH - V.*((pars.M.*pars.phi)'*Ndeb') ;
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
odefun = @(t, y) seivd(t, y, pars,id);

% Run the ODE solver
[t_out, y_out] = ode45(odefun, tvec, y0);



%% total bacteria
S = y_out(:,id.S);

for i = 1:pars.NH 
    B_out(:,i) = S(:,i) + sum(y_out(:,id.Emat(i,:)),2) + sum(y_out(:,id.Imat(i,:)),2);
end


%% virus

V_out = y_out(:,id.V);

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



subplot(2,5,4)
errorbar(time/60,mean(1e3*host4'),std(1e3*host4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20)
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e9]);
    xlim([0 16]);
  xticks([0 2 4 6 8 10 12 14 16]);
  axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA H100','FontSize',18);
    
    plot(t_out,B_out(:,4),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,5)
errorbar(time/60,mean(1e3*host5'),std(1e3*host5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
ylim([1e5 1e9]);
    xlim([0 16]);
    xticks([0 2 4 6 8 10 12 14 16]);
    axis('square');
    yticks([1e5 1e6 1e7 1e8]);
    title('PSA 13-15','FontSize',18);
    
    plot(t_out,B_out(:,5),'-','Color',color_ofthe_fit,'LineWidth',linewidth);

%xlabel("Time (hours)");
%ylabel("Host density (cell/ml)");






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

subplot(2,5,9)
errorbar(time/60,mean(1e3*virus4'),std(1e3*virus4'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
   xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
   yticks([1e4 1e6 1e8 1e10]);
   title('PSA HP1','FontSize',18);
    
    plot(t_out,V_out(:,4),'-','Color',color_ofthe_fit,'LineWidth',linewidth);




subplot(2,5,10)
errorbar(time/60,mean(1e3*virus5'),std(1e3*virus5'),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[70/255,130/255,180/255],Color=[70/255,130/255,180/255]);hold on;
set(gca, 'YScale', 'log');set(gca,'fontname','times')  % Set it to times
ylim([1e4 1e11]);
    xlim([0 16]);
 xticks([0 2 4 6 8 10 12 14 16]);
    set(gca,'FontSize',20);
    axis('square');
  yticks([1e4 1e6 1e8 1e10]);
  title('PSA HS6','FontSize',18);
    
    plot(t_out,V_out(:,5),'-','Color',color_ofthe_fit,'LineWidth',linewidth);
    %legend('Data','95% confidence interval','Bayesian fit');
    %legend('Box','off');

han=axes(hf4,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
set(gca,'FontSize',20);
set(gca,'fontname','times')  % Set it to times
xlabel("Time (hours)");

