% @ CopyRight Mohammed Alamsar
clc,clear,close all


% dinfo = dir('100ms/*.mat');               
dinfo = dir('10sec/*.mat');    %%%%%%%%%%%%%%%%%%%%%%%%%%%             

thisfilename=cell(1,length(dinfo));
for K = 1 : length(dinfo)
    thisfilename{K} = dinfo(K).name;  %just the name
end


for j=0:23
   traces{j+1} = [num2str(j),'.mat']; 
end


 ResultsDirectory = [ '10sec/'];  %%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(traces)
    traceIndex=i
    traceName  =  traces{i};
    
    traceNameN = ['data_rate_',num2str(i-1)]; %  erase(traceName,'.mat') 
    inputMatFile = ['10sec/' traceName];             %%%%%%%%%%%%%%%%%%%%%%%%%%%

    stru=load(inputMatFile);
    dataRate= getfield(stru, traceNameN);
    
    
    %% STATISTICS Values
    variance =var(dataRate(1:end-1))  ; % bit^2
    ave = mean(dataRate);
    sigma = std(dataRate);
    aggDataLenght= length(dataRate);
    maxValue= max(dataRate);
    
    
    %% Step 1 estimate  x_min and alpha
    [alpha, xmin, ntail] = plfit(dataRate );
    
    %% Step 2 estimates the uncertainty in the estimated power-law parameters
     %alphaUn =0 ;  xminUn=0; ntailUn= 0; 

    [alphaUn, xminUn, ntailUn]= plvar(dataRate,'sample',100);
    %p =plpva(x, xmin )
 
    %% store the results in a csv file
    c = {'traceId' , 'n' , 'xmax' , 'mean', 'std' , 'alpha' ,'xmin' , 'ntail' ,'alphaUn' , 'xminUn' ,'ntailUn' ;[traceIndex] , [aggDataLenght], [maxValue], [ave] ,[sigma], [alpha], [xmin],[ntail], [alphaUn], [xminUn],[ntailUn] };
    
    CSVName=[ResultsDirectory, traceNameN, 'Results' ,'.csv']   ;
    
    fid = fopen(CSVName,'w');
    fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',c{1,:});
    fprintf(fid,'%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f ',c{2:end,:});
    fclose(fid);
    
    
    % CDF plot  (on log axes) data rates vs  power-law distribution
    figure, plplot(dataRate, xmin, alpha);
    figlogCDFName=[ResultsDirectory, traceNameN, 'logCDF' ,'.fig']
    savefig(figlogCDFName);
end
