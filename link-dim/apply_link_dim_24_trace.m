% @ CopyRight Mohammed Alamsar
clc,clear,close all

epsi =0.5  ;
traceLoc=['1sec/'];
T=1;


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


for i=1:length(traces)
    traceIndex=i
    traceName  =  traces{i};
    
    traceNameN = erase(traceName,'.mat')
    inputMatFile = [traceLoc  traceName];
    stru=load(inputMatFile);
    
    dataRateVar = ['data_rate_',num2str(i-1)]
    dataRate= getfield(stru, dataRateVar); % bits/sec
    
    
    aggData= dataRate*T;  % bit
    
    %% STATISTICS Values
    VAR=var(aggData(1:end-1))  ; % bit^2
    ave=mean(aggData(1:end-1))/T ; %% bit/sec
    sigma =std(aggData(1:end-1));  %% bit
    aggDataLenght= length(aggData);
    maxValue= max(dataRate);
     

    %% Finding the Value of the Capacity
    %%% first equ --- Remco equ
    C1=ave+((1/T)*sqrt(-2*log(epsi)*VAR));
    % safety_margin(kk)=(1/T)*sqrt(-2*log(epsi(kk))*VAR(kk));
    % safety_margin_part(kk)=(1/T)*sqrt(VAR(kk));
        
    %%% second equ --- Nick equ
    C2=fsolve('solveC2K',1.4*C1);
    
    %%% third equ --- direct Gaussian way
    C3=icdf('norm',1-epsi,ave,sigma) ;
    
    %%% fourth equ --- Generalized Extreme Value
    GEVfitRes=gevfit(dataRate);
    C4=icdf('Generalized Extreme Value',1-epsi,GEVfitRes(1),GEVfitRes(2),GEVfitRes(3));
    
    %%% Fifth equ --- Log-normal
    LOGnormalfit=lognfit(dataRate);
    C5=icdf('Lognormal',1-epsi,LOGnormalfit(1),LOGnormalfit(2));
    
    epsi_c1(i)=length(find(dataRate(1:end-1)>C1))/(length(dataRate)-1);  % remco
    epsi_c2(i)=length(find(dataRate(1:end-1)>C2))/(length(dataRate)-1);  % nick
    epsi_c3(i)=length(find(dataRate(1:end-1)>C3))/(length(dataRate)-1);  % direct
    epsi_c4(i)=length(find(dataRate(1:end-1)>C4))/(length(dataRate)-1);  % GEV
    epsi_c5(i)=length(find(dataRate(1:end-1)>C5))/(length(dataRate)-1);  % lognormal
    
end



%% avg(epsi)
avgEpsiC1 = mean(epsi_c1)
avgEpsiC2 = mean(epsi_c2)
avgEpsiC3 = mean(epsi_c3)
avgEpsiC4 = mean(epsi_c4)
avgEpsiC5 = mean(epsi_c5)

%% avg(epsi-emp)
avgEpsiEmpC1 = mean(abs(epsi_c1-epsi))
avgEpsiEmpC2 = mean(abs(epsi_c2-epsi))
avgEpsiEmpC3 = mean(abs(epsi_c3-epsi))
avgEpsiEmpC4 = mean(abs(epsi_c4-epsi))
avgEpsiEmpC5 = mean(abs(epsi_c5-epsi))

%% stderr(epsi-emp)
stderrEpsiEmpC1 = std(abs(epsi_c1-epsi))
stderrEpsiEmpC2 = std(abs(epsi_c2-epsi))
stderrEpsiEmpC3 = std(abs(epsi_c3-epsi))
stderrEpsiEmpC4 = std(abs(epsi_c4-epsi))
stderrEpsiEmpC5 = std(abs(epsi_c5-epsi))


res= [
    T, epsi, avgEpsiC1 , avgEpsiEmpC1 , stderrEpsiEmpC1;
    T, epsi, avgEpsiC2 , avgEpsiEmpC2 , stderrEpsiEmpC2;
    T, epsi, avgEpsiC3 , avgEpsiEmpC3 , stderrEpsiEmpC3;
    T, epsi, avgEpsiC4 , avgEpsiEmpC4 , stderrEpsiEmpC4;
    T, epsi, avgEpsiC5 , avgEpsiEmpC5 , stderrEpsiEmpC5
    ]
