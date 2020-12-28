% @Mohammed Alasmar 
clc,clear, close all

MarkerSizeA=14, LineWidthA=1;

x =  [1:24]; %  24 hour 
load('LR24MawiPythonResults100msec.mat')

lrlog = mawi24HourLr100msecLognorm;
lrexp = mawi24HourLr100msecExp;
lrwei = mawi24HourLr100msecWeibull;

plog= mawi24HourPvalue100msecLognorm;
pexp = mawi24HourPvalue100msecExp; 
pwei = mawi24HourPvalue100msecWeibull; 
 

plot(x,zeros(1,length(x)), '-g')

%% Weibull
for i = 1:length(x)
    [lrWeSort(i),ix] = min(lrwei);       % Find Minimum & Index
    lrwei(ix) = [1000];                  % Set Previous Minimum To Empty
    pWeiSort(i)=pwei(ix)
end

hold on
h1 =plot(x,lrWeSort,'--d',   'Color', 'm','MarkerSize', MarkerSizeA , 'LineWidth',LineWidthA)
h1 =plot(x,lrWeSort,'--d',   'Color', 'm','MarkerSize', MarkerSizeA , 'LineWidth',LineWidthA)

hold on

for i=1:length(x)
      if pwei(i)>0.1 
          plot(i,lrWeSort(i), 'om', 'MarkerSize', MarkerSizeA);
      end 
end

%% lognormal
for i = 1:length(x)
    [lrlogSort(i),ix] = min(lrlog);       % Find Minimum & Index
    lrlog(ix) = [1000];                   % Set Previous Minimum To Empty
    plogSort(i)=plog(ix)
end

hold on
h2 =plot(x,lrlogSort,'--x',   'Color', 'k','MarkerSize', MarkerSizeA , 'LineWidth',LineWidthA)
hold on

for i=1:length(x)
      if plogSort(i)>0.1 
          plot(i,lrlogSort(i), 'ok', 'MarkerSize', MarkerSizeA);
      end 
end


%% Exp
for i = 1:length(x)
    [lrExpSort(i),ix] = min(lrexp);       % Find Minimum & Index
    lrexp(ix) = [1000];                   % Set Previous Minimum To Empty
    pExpSort(i)=pexp(ix)
end

hold on
h3 =plot(x,lrExpSort,'--+',   'Color', 'b','MarkerSize', MarkerSizeA , 'LineWidth',LineWidthA)
hold on

for i=1:length(x)
      if pExpSort(i)>0.1 
          plot(i,lrExpSort(i), 'ob', 'MarkerSize', MarkerSizeA);
      end 
end


[lgd, icons, plots, txt] =  legend([h1 h2 h3],{   'Weibull' , 'Log-normal', 'Exponential'},'FontSize',26)
%title(tit)  %%%%
xlabel('Rank of trace id (1-hour long each)')
ylabel('Normalised  LLR')
% xlim([60 107.5])
grid on ,
set(gca,'fontsize',26)
set(gca,'FontName','Times')
savefig('in.fig')
h3 =openfig('in.fig')
 
xlim([1 24])
set(gca, 'XTick', 0:4:24);

grid on ,
set(gca,'fontsize',28)
set(gca,'FontName','Times')
savefig('in.fig')
h3 =openfig('in.fig')
box on;
set(h3,'Units','Inches');
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
 
