clc,clear ,close all


% load('mawiPL100msecresults.mat')
%  x =  [1:110];
% tit=['107 Mawi traces'] 
% lrlog = mawiLr100msecLognorm;     lrlog(108:110) = lrlog(100:102);
% lrexp = mawiLr100msecExp;        lrexp(108:110) = lrexp(100:102);
% lrwei = mawiLr100msecWeibull;    lrwei(108:110) = lrwei(100:102);
% plog= mawiPvalue100msecLognorm;  plog(108:110) = plog(100:102);
% pexp = mawiPvalue100msecExp;      pexp(108:110) = pexp(100:102);  
% pwei = mawiPvalue100msecWeibull;   pwei(108:110) = pwei(100:102);

% con= 5;
 

% x =  [1:40];
% tit=['40 Twente traces']
%   load('twentePL100msecresults.mat')
% lrlog = twenteLr100msecLognorm;
% lrexp = twenteLr100msecExp;
% lrwei = twenteLr100msecWeibull;
% plog= twentePvalue100msecLognorm;
% pexp = twentePvalue100msecExp;
% pwei = twentePvalue100msecWeibull;
% con= 4;
%  
% x =  [1:25];
% tit=['25 Auckland traces']
% load('aucklandPL100msecresults.mat')
% lrlog = aucklandLr100msecLognorm;
% lrexp = aucklandLr100msecExp;
% lrwei = aucklandLr100msecWeibull;
% plog= aucklandPvalue100msecLognorm;
% pexp = aucklandPvalue100msecExp;
% pwei = aucklandPvalue100msecWeibull;
% con= 3;


% 
% x =  [1:30];
% tit=['30 Waikato traces'] 
% load('waikatoPL10msecResults.mat')
% lrlog = waikatoLr10msecLognorm;
% lrexp = waikatoLr10msecExp;
% lrwei = waikatoLr10msecWeibull;
% plog= waikatoPvalue10msecLognorm;
% pexp =waikatoPvalue10msecExp; 
% pwei =waikatoPvalue10msecWeibull; 
% con= 2;


x =  [1:27];
tit=['27 Caida traces']
load('caidaPL10msecResults.mat')
lrlog = ispdslLr10msecLognorm;
lrexp = ispdslLr10msecExp;
lrwei = ispdslLr10msecWeibull;
plog= ispdslPvalue10msecLognorm;
pexp = ispdslPvalue10msecExp; 
pwei = ispdslPvalue10msecWeibull; 
con= 1;





for i = 1:length(x)
    [lrlogSort(i),ix] = min(lrlog);       % Find Minimum & Index
    lrlog(ix) = [1000];                 % Set Previous Minimum To Empty
    plogSort(i)=plog(ix)
end
 
 
 
for i = 1:length(x)
    [lrexpSort(i),ix] = min(lrexp);       % Find Minimum & Index
    lrexp(ix) = [1000];                 % Set Previous Minimum To Empty
    pexpSort(i)=pexp(ix)
end

for i = 1:length(x)
    [lrweiSort(i),ix] = min(lrwei);       % Find Minimum & Index
    lrwei(ix) = [1000];                 % Set Previous Minimum To Empty
    pweiSort(i)=pwei(ix)
end

 


 plot(x,zeros(1,length(x)), '-g')
hold on
%errorbar(x,lrweiSort,pweiSort,'--',  'LineStyle', '-', 'Color', 'k')
h1 =  plot(x,lrweiSort, '--d', 'Color', 'm','MarkerSize', aa , 'LineWidth',bb)
%errorbar(x,lrlogSort,plogSort,'.',  'LineStyle', '-', 'Color', 'r')
hold on
h2 =plot(x,lrlogSort,'--x',   'Color', 'k','MarkerSize', aa , 'LineWidth',bb)

hold on
%errorbar(x,lrexpSort,pexpSort,'*',  'LineStyle', '-', 'Color', 'b')
h3 =plot(x,lrexpSort,'--+', 'Color', 'b','MarkerSize', aa , 'LineWidth',bb)

 



for i=1:length(x)
      if pexpSort(i)>0.1 
          plot(i,lrexpSort(i), 'ob', 'MarkerSize', cc);
      end 
end


for i=1:length(x)
      if pweiSort(i)>0.1 
          plot(i,lrweiSort(i), 'om', 'MarkerSize', cc);
      end 
end

 for i=1:length(x)
      if plogSort(i)>0.1 
          plot(i,lrlogSort(i), 'ok', 'MarkerSize', cc2);
      end 
 end

[lgd, icons, plots, txt] =  legend([h1 h2 h3],{   'Weibull' , 'Lognormal', 'Exponential'},'FontSize',26)
%title(tit)  %%%%
xlabel('Rank of trace')
ylabel('Normalised   LLR')

if con == 1
xlim([1 27.2])
end

if con == 2
xlim([1 30.2])
end

if con == 3
xlim([1 25.2])
ylim([-50 10])

end

if con == 4
xlim([1 40.2])
  ylim([-70 20])
end


if con == 5
 xlim([60 110])
  ylim([-7 10])
end
% xlim([60 107.5])
grid on ,
set(gca,'fontsize',28)
set(gca,'FontName','Times')
savefig('in.fig')
h3 =openfig('in.fig')

set(h3,'Units','Inches');
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%saveas(h3,'incast.pdf')

