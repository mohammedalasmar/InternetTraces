%   @Mohammed Alasmar
%       clc,clear,close all
close all
%  load('nintyfifthCaida.mat')
%  load('ninetyfifthWaikato.mat')
%  load('ninetyfifthAuckland.mat')
   load('nintyfiftLong.mat')
%  load('ninetyfifthTwente.mat')


close all

den= 10^6;
Clol=C5/den;    % log
Cnorm=C3/den;   % actual norm
Cgev=C4/den;    % gev
Cheavy=C7/den;  % pareto

%%% caida
xlim([2 5.6])
ylim([2 6.4])
%%% Waikato
xlim([0 110])
ylim([0 130])
%%% Auckland
xlim([12 165])
ylim([0 252])
%%% Twente
xlim([0 35])
ylim([0 52])

%%% Mawi
xlim([100 850])
ylim([110 1350])


Cmeen= C1/den;
perAct = per/den;
aa=1.3;
bb= 12;
 


plot(perAct,perAct,'r','LineWidth',2)
hold on, h3= plot(perAct,Cmeen, '*b', 'MarkerSize',bb,'LineWidth',aa)
hold on, h2= plot(perAct,Cheavy,'.g','colo', [0, 0.74, 0],'MarkerSize',4*bb,'LineWidth',aa)
hold on, h1= plot(perAct,Cgev, 'ok','MarkerSize',bb,'LineWidth',1.2*aa)
plot(perAct,perAct,'r','LineWidth',2)


ylabel('Predicted value (Mbps)')
xlabel('Actual value (Mbps)')
[lgd, icons, plots, txt] =  legend([h1 h2 h3 ], {'Log-normal', 'Weibull', 'Gaussian'},'FontSize',30 )


grid on
set(gca,'fontsize',32)
set(gca,'FontName','Times')
savefig('in.fig')
h3 =openfig('in.fig')
set(h3,'Units','Inches');
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(h3,'test.pdf')
