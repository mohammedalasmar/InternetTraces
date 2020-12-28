
clc,clear, close all


traceLoc=['100ms/']; aa=24;

T = 0.1;
traceLocMat = [traceLoc '*.mat']
dinfo = dir(traceLocMat);

thisfilename=cell(1,length(dinfo));
for K = 1 : length(dinfo)
    thisfilename{K} = dinfo(K).name;  %just the name
end

traces = thisfilename;

for j=0:23
   traces{j+1} = [num2str(j),'.mat']; 
end


for i=1:aa
    traceIndex=i
    traceName  =  traces{i};

    traceNameN = ['data_rate_',num2str(i-1)]; %  erase(traceName,'.mat') 
    inputMatFile = [traceLoc  traceName];
    stru=load(inputMatFile);
    dataRate= getfield(stru, traceNameN);
    aggData= dataRate*T;  % bit
    aggregationTime =  T;
    %% STATISTICS Values
    variance =var(dataRate(1:end-1))  ; % bit^2
    ave = mean(dataRate);
    sigma = std(dataRate);
    aggDataLenght= length(dataRate);
    maxValue= max(dataRate);
    statistics={'aggDataLenght', num2str(aggDataLenght) ;'mean', num2str(ave) ;'std', num2str(sigma) ;'max', num2str(maxValue) };
    x = dataRate;
 
    %%
    dataRateSort= sort(dataRate);
    %tNew= linspace(0,100, length(t))
    %figure, bar(tNew, dataRateSort, 'k')
    xlabel('percentile'), ylabel('data rate (bps)')
    xlim([0 100])
    ylim([min(dataRate)-1000000 max(dataRate)+1000000 ])
    
    yp=[1:1000:max(dataRate)];
    xp=95*ones(1,length(yp));
    % hold on, plot(xp,yp ,'-r')
    
    per(i)= prctile(dataRate,95)
    
    C1(i)=ave+(1/aggregationTime)*sqrt(-2*log(1-0.95)*variance);
    C3(i)=icdf('norm',0.95,ave,sigma/T);
    GEVfitRes=gevfit(dataRate);
    C4(i)=icdf('Generalized Extreme Value',0.95,GEVfitRes(1),GEVfitRes(2),GEVfitRes(3));
    LOGnormalfit=lognfit(dataRate);
    C5(i)=icdf('Lognormal',0.95,LOGnormalfit(1),LOGnormalfit(2));
    weibullfit=wblfit(dataRate);
    C6(i)=icdf('Weibull',0.95,weibullfit(1),weibullfit(2));
    
    paretofit=gpfit(dataRate);
    C7(i)=icdf('Generalized Pareto',0.95,paretofit(1),paretofit(2));
end


 
