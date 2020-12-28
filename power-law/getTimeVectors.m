
clc,clear,close all

% load all the txt files
dinfo = dir('all/*.txt');
thisfilename=cell(1,length(dinfo));
for K = 1 : length(dinfo)
    thisfilename{K} = dinfo(K).name;  %just the name
end

traces = thisfilename;

traceLoc=['all/'];



for traceID=1:length(traces)
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
    
    clear Struc1;
    timeVec=['timeVecTrace-', int2str(traceID) ,'.mat']
    newName = ['timeVecTrace' int2str(traceID)];
    Struc1.(newName) = timeValues;
    save(timeVec, '-struct', 'Struc1')

end
