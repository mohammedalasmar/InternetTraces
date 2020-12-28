% @Mohammed Alasmar
clc,clear,close all
load('mawiTrace16dataRateAtAggTime100msec.mat')

T = 0.1;
aggregationTimeScale = 0.1;
%  getting data rates from the input trace
aggregationTime = aggregationTimeScale; 
dataRate= mawiTrace16dataRateAtAggTime100; 

 
%  plot data rates PDF
figure , histogram(dataRate/10^9,300,'Normalization','probability', 'FaceColor' , 'k')
xlabel('Data rate (Gbps)'), 
grid on
legend('PDF')
legend boxoff
set(gca,'fontsize',45)
set(gca,'FontName','Times')
savefig('in.fig')
h3 =openfig('in.fig')

set(h3,'Units','Inches');
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])



load('mawiTrace12dataRateAtAggTime100msec.mat') % ?anomalous
aggregationTimeScale=0.1;
T = aggregationTimeScale;
%  getting data rates from the input trace
aggregationTime = aggregationTimeScale; 
dataRate= mawiTrace12dataRateAtAggTime100; 

 
 

%  plot data rates PDF
figure , histogram(dataRate/10^9,300,'Normalization','probability', 'FaceColor' , 'k')
xlabel('Data rate (Gbps)'), 
%ylabel('PDF')
legend('PDF')
legend boxoff

 grid on
%ax = gca; ax.XAxis.Exponent = 6
set(gca,'fontsize',45)
set(gca,'FontName','Times')
savefig('in.fig')
h3 =openfig('in.fig')

set(h3,'Units','Inches');
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

 
 


