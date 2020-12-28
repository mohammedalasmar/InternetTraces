clc,clear,close all

load('rate_0.1_sec.mat')
load('rate_0.5_sec.mat')
load('rate_1_sec.mat')
load('rate_10_sec.mat')
load('rate_20_sec.mat')
load('rate_60_sec.mat')

 
t24 = t24_0_1_sec;  all_data_rate = all_data_rate_0_1_sec; plot(t24/(3600),all_data_rate/10^9, 'r', 'LineWidth', 2)
% hold on, t24 = t24_0_5_sec;  all_data_rate = all_data_rate_0_5_sec; plot(t24/(3600),all_data_rate/10^9, 'c')
hold on, t24 = t24_1_sec;  all_data_rate = all_data_rate_1_sec; plot(t24/(3600),all_data_rate/10^9, 'g', 'LineWidth', 2)
% hold on, t24 = t24_10_sec;  all_data_rate = all_data_rate_10_sec; plot(t24/(3600),all_data_rate/10^9, 'k')
% hold on, t24 = t24_20_sec;  all_data_rate = all_data_rate_20_sec; plot(t24/(3600),all_data_rate/10^9, 'y')
hold on, t24 = t24_60_sec;  all_data_rate = all_data_rate_60_sec; plot(t24/(3600),all_data_rate/10^9, 'b', 'LineWidth', 2)

 lgd =  legend({ 'T=100ms ','T= 1s ','T=60s '},'FontSize',35 , 'Orientation' , 'horizontal', 'NumColumns',3)
ylim([0 2.2])
xlim([0 24])
set(gca,'XTick',[0:1:24]) 
xlabel('Time (hour)')
ylabel('Data rate (Gbps)')
% xticks(0:2:24)
set(gca,'fontsize',35)
set(gca,'FontName','Times')
savefig('in.fig')
h3 =openfig('in.fig')
set(h3,'Units','Inches');
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
box on;
%  autocorr(all_data_rate_0_1_sec(100000:200000),'NumLags',2022)






figure,
t24 = t24_1_sec;  all_data_rate = all_data_rate_1_sec; histogram(all_data_rate/10^9,300,'FaceColor','k','Normalization','probability')
% hold on, t24 = t24_60_sec;  all_data_rate = all_data_rate_60_sec;histogram(all_data_rate/10^9,80,'FaceColor','g','Normalization','probability')
grid on
 lgd =  legend({ 'T= 1s'},'FontSize',25 , 'Orientation' , 'horizontal', 'NumColumns',3)
% ylim([0 2.2])
xlabel('Data rate (Gbps)')
ylabel('PDF')
% xticks(0:2:24)
set(gca,'fontsize',30)
set(gca,'FontName','Times')
savefig('in.fig')
h3 =openfig('in.fig')
set(h3,'Units','Inches');
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
box on;






figure,
t24 = t24_0_1_sec;  all_data_rate = all_data_rate_0_1_sec;
h1=cdfplot(all_data_rate/10^9), set( h1, 'LineStyle', '-', 'Color', 'b');
hold on
t24 = t24_60_sec;  all_data_rate = all_data_rate_60_sec;
h1=cdfplot(all_data_rate/10^9), set( h1, 'LineStyle', '-', 'Color', 'r');



 %%%%

dataRate = all_data_rate_1_sec;

lenn = length(dataRate);
lags_adf = floor(12*(lenn/100)^(1/4));
lags_kspp = floor(lenn^(1/2));
lags_pp = floor(lenn^(1/4));
lags_kspp = floor(lenn^(0.6));

diff_data = diff(dataRate,1);
 
[h_kpss,pValue_kpss , stat_kpss,cValue_kpss] = kpsstest(dataRate, 'trend' , false , 'alpha',0.05, 'lags' , lags_kspp ); % 'alpha' significance levels for the hypothesis tests. 95% confidence interval
[h_kpss_trend,pValue_kpss_trend , stat_kpss_trend,cValue_kpss_trend] = kpsstest(dataRate, 'trend' , true , 'alpha',0.05, 'lags' , lags_kspp); % 'alpha' significance levels for the hypothesis tests. 95% confidence interval
[h_kpss_diff,pValue_kpss_diff , stat_kpss_diff,cValue_kpss_diff] = kpsstest(diff_data, 'trend' , false , 'alpha',0.05, 'lags' , lags_kspp); % 'alpha' significance levels for the hypothesis tests. 95% confidence interval

% [h_adf,pValue_adf , stat_adf,cValue_adf] = adftest(dataRate, 'alpha',0.05 ,'model',{'AR'}, 'lags' ,lags_adf-13);
% [h_pp,pValue_pp , stat_pp ,cValue_pp ] = pptest(dataRate, 'alpha',0.05, 'model',{'ARD'},'lags' , lags_pp);


[h_adf,pValue_adf , stat_adf,cValue_adf] = adftest(dataRate(1:5:end), 'alpha',0.05 ,'model',{'AR'}, 'lags' , 40) 
[h_pp,pValue_pp , stat_pp ,cValue_pp ] = pptest(dataRate(1:200:end), 'alpha',0.05, 'model',{'AR'},'lags' ,7) 


Kpss_combine = max([pValue_kpss,pValue_kpss_trend])
all_p_values = [pValue_adf, pValue_pp, Kpss_combine]
all_p_values_1_sec = all_p_values
%  save('all_p_values_24_hr.mat','all_p_values_0_1_sec','all_p_values_0_5_sec','all_p_values_1_sec','all_p_values_10_sec','all_p_values_20_sec','all_p_values_60_sec')




%%$

load('all_p_values_24_hr.mat')
all_adf_pvalues = [all_p_values_0_1_sec(1), all_p_values_0_5_sec(1), all_p_values_1_sec(1), all_p_values_10_sec(1), all_p_values_20_sec(1), all_p_values_60_sec(1)]
all_pp_pvalues = [all_p_values_0_1_sec(2), all_p_values_0_5_sec(2), all_p_values_1_sec(2), all_p_values_10_sec(2), all_p_values_20_sec(2), all_p_values_60_sec(2)]
all_kpss_pvalues = [all_p_values_0_1_sec(3), all_p_values_0_5_sec(3), all_p_values_1_sec(3), all_p_values_10_sec(3), all_p_values_20_sec(3), all_p_values_60_sec(3)]



figure

r=[ 51 102  0 
    102 204 0
    128 255 0
    178 255 102
    229 255 204
    244 244 244
   ];
r = r/255;
 
z = [all_adf_pvalues;all_pp_pvalues;all_kpss_pvalues];
errBar = zeros(size(z))
errorbar_groups(z',errBar','bar_colors' , r ,'bar_width',0.8,'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',1.5}, 'errorbar_width',3, 'bar_names' , {'ADF' ,'PP', 'KPSS'} )
% [lgd, icons, plots, txt] =  legend({ 'T=100ms ','T= 500ms ','T=1s ',  'T=10s ' , 'T=20s ', 'T=60s '},'FontSize',22 , 'Orientation' , 'horizontal', 'NumColumns',3)
% [lgd, icons, plots, txt] =  legend({ 'T=100ms ','T= 500ms ','T=1s '},'FontSize',20 , 'Orientation' , 'horizontal', 'NumColumns',3)
 lgd =  legend({ 'T=100ms ','T= 500ms ','T=1s ',  'T=10s ' , 'T=20s ', 'T=60s '},'FontSize',19 , 'Orientation' , 'horizontal', 'NumColumns',3)
% legend boxoff
% lgd = legend;
%   lgd.NumColumns = 2;
% lgd.Orientation = 'vertical';

grid on
xlabel('Stationarity tests')
ylabel('p-value')
% xticks(0:2:24)
set(gca,'fontsize',30)
set(gca,'FontName','Times')
savefig('in.fig')
h3 =openfig('in.fig')
set(h3,'Units','Inches');
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
box on;
% x = [1:3];
% y = [1:6]; 
% z = [all_adf_pvalues;all_pp_pvalues;all_kpss_pvalues];
% [x, y] = meshgrid(1:3,1:6);
% figure
% hh= imagesc(flipud(z')), 
% % axis equal tight, 
% colorbar
% h = colorbar;
%  colormap(gray)
%    colorMap = jet(444);
% %  colormap(gray)
% 
% %  colormap(hsv(122))
% %  colormap(colorMap)
% 
% set(gca, 'YTick', 1:9, 'YTickLabel', 9:-1:1);
% box on; set(gca, 'FontSize',14) 
% names = {'60s';'20s';'10s'; '1s' ; '500ms';'100ms'};  set(gca,'ytick',[1:6],'yticklabel',names)
% xlabel('ADF test')
% ylabel('Timescale')
% title('p-value')



