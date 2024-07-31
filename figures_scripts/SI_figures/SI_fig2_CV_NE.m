clear all;
clc;
addpath(genpath(pwd));

mean_tau = 2;
x = 0:0.01:2*mean_tau;

i=0;
NE = [1,5,10,25,100];


hf = figure;

for shape = NE
rate = shape/mean_tau;
pdf = ((rate^shape) * (x.^(shape -1)) .* exp(-rate*x)/gamma(shape));
i=i+1;
ha(i) = subplot(length(NE),2,2*i-1)
plot(x,pdf,'LineStyle','-','LineWidth',2);hold on;yticks([]);
set(gca,'FontSize',20,'FontName','Times');

line([mean_tau, mean_tau], [0, max(pdf)], 'Color', [0.1,0.1,0.1],'LineStyle','--',LineWidth=2);

if i == 5
    xlabel('\tau (hrs)');
elseif i == 3   
    ylabel('PDF of \tau');
end


CV = 1/sqrt(shape);
str = strcat('cv = ',num2str(CV),' N_E =', num2str(shape));
yticks([0 0.3 0.6 0.9 1.2 1.5 1.8 2.1]);
%annotation('textbox',dim,'String',str,'FitBoxToText','on');

%ha(i).Position = ha(i).Position + [0 0 0 0.5]; 

end

set(gca,'fontname','times')



pos = get(ha, 'position');
dim = cellfun(@(x) x.*[1 1 0.5 0.5], pos, 'uni',0);

% 2 points after decimal


for i = 1: length(NE)

CV = 1/sqrt(NE(i));
annotation(hf, 'textbox', dim{i}, 'String',  strcat('CV = ','  ',num2str(CV,2), '; ',' N_E =   ',num2str(NE(i)) )  ,'Position', ha(i).Position + [0.2 0 0 0],'EdgeColor','none', 'FontSize',20,'verticalalignment', 'top','FitBoxToText','on','FontName','times') ;

end

%% functions inclusion
time_free_phages = 0:0.01:4;

for i = 1:length(NE)

NE_cases = NE(i); %number of exposed class

moi_mean = 0.1;

S0 = 1e8;
V0 = S0*moi_mean;

y(1,1) = S0;
y(1,2:NE_cases+1) = 0;
y(1,NE_cases+2) = 0 ;
y(1,NE_cases+3) = V0;

theta_raunak = [0.2,13e-8,2,200];
dilution_factor = 100;
[time2,y_series2,time_abs,pre_dil] = one_step_simulate_particular_points(time_free_phages,y,theta_raunak,NE_cases, dilution_factor);

subplot(length(NE),2,2*i);
semilogy(time2,y_series2(end,:),'-k', 'LineWidth',2);

if i == 5
    xlabel('Time (hrs)');
elseif i == 3  
    ylabel('Phage density (virions/ml)');
end
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10]);

set(gca,'FontSize',20)
set(gca,'fontname','times')
end


%%

han=axes(hf,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
%ylabel(han,'PDF of \tau ');
%xlabel(han,'\tau (hours)');
set(gca,'FontSize',20)

set(gca,'fontname','times')
