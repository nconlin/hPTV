function CC = TrackDatatoCC(tracks)

maxTime = -1;
minTime = -1;
for ii = 1:numel(tracks)
    newMin = min(tracks(ii).time);
    newMax = max(tracks(ii).time);
    if newMin < minTime
        minTime = newMin;
    end
    if newMax > maxTime
        maxTime = newMax;
    end
end
blnkCC = struct('X',[],'Y',[]);
CC = repmat(blnkCC,[maxTime-minTime+1 1]);

for ii = 1:numel(tracks)
    times = tracks(ii).time;
    timeCount = 1;
    for jj = times
        % append current positions
        CC(jj).X = [CC(jj).X tracks(ii).pos(timeCount,1)]; 
        CC(jj).Y = [CC(jj).Y tracks(ii).pos(timeCount,2)];
        timeCount = timeCount + 1;
    end
end
