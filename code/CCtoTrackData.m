function [X,T] = CCtoTrackData(CC)

X = [];
T = [];
nFrames = numel(CC);
for ii = 1:nFrames
    xnow = CC(ii).X;
    ynow = CC(ii).Y;
    nParts = length(xnow);
    if size(xnow,1) == 1
        xnow = xnow';
        ynow = ynow';
    end
    X = [X; xnow ynow];
    T = [T; ii*ones(nParts,1)];
end
