% @ CopyRight Mohammed Alamsar

clc,clear,close all

% load all the txt files
dinfo = dir('all/*.txt');
thisfilename=cell(1,length(dinfo));
for K = 1 : length(dinfo)
    thisfilename{K} = dinfo(K).name;  %just the name
end

traces = thisfilename;

traceLoc=['all/'];

 aggregationTimeScale = [0.005 0.01 0.025 0.05 0.1 0.5 1 2 5];
%aggregationTimeScale = [0.01];


for traceID=1:length(traces)

for timeIter = 1:length(aggregationTimeScale)

    traceIndex=traceID
    traceName  =  traces{traceID};
    
    
    
    %% get packets lengths and their times from the trace
    %traceLoc=['../traces/mawiTraces/group1-19/' traceName ]
    trace=[traceLoc traceName]
    fid = fopen(trace);
    timePktLenMatrix = fscanf(fid,'%g %g',[2 inf]);
    timePktLenMatrix = timePktLenMatrix';
    fclose(fid);
    timeValues=timePktLenMatrix(:,1); % time vector
    pktsLength=timePktLenMatrix(:,2)*8; % bits
    
    traceNameN = erase(traceName,'.txt');
    
    
    %% getting data rates from the input trace
    aggregationTime = aggregationTimeScale(timeIter) 
    pickedTimesIndex=0;
    dataRate = 0;
    aggData = 0;
    j=1;
    T = aggregationTime;
    
    % aggregation based on selected aggregationTime
    for i=1:length(timeValues)
        if timeValues(i)>=T
            pickedTimesIndex(j)=i;
            T = T + aggregationTime;
            j=j+1;
        end
    end
    
    pickedTimesIndex=[1 pickedTimesIndex length(timeValues)+1];
    
    % aggreate packets
    for u=1:length(pickedTimesIndex)-1
        aggData(u)=sum(pktsLength(pickedTimesIndex(u):pickedTimesIndex(u+1)-1));
    end
    
    % data rate
    dataRate=aggData/aggregationTime; % bits/sec


    clear Struc1;
    timeScaleMs = aggregationTimeScale(timeIter)*1000;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dataRateAtAggTime=['twenteTrace' , int2str(traceIndex), 'dataRateAtAggTime', int2str(timeScaleMs) ,'msec.mat']
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newName = ['twenteTrace' int2str(traceIndex)  'dataRateAtAggTime' int2str(timeScaleMs)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Struc1.(newName) = dataRate;
    save(dataRateAtAggTime, '-struct', 'Struc1')
end
end
