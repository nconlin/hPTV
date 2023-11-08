function [X,T,E] = readSTM(file,frames)
% read data from running c version of stm 

if numel(frames) == 1
    frames = 1:frames(1);
elseif numel(frames) == 2
    frames = frames(1):frames(2);
end
nframes = length(frames);

% allocate some space
X = []; T = []; E = [];
for ii = 1:nframes

    % read in data for this frame
    dataPath = ['/frame' num2str(frames(ii)) '_xyze'];
    data = h5read(file,dataPath);
    
    % append x,y,z positions
    X = [X; data(:,1:3)];

    % append error
    E = [E; data(:,4)];

    % set the time as the same for each particle this frame
    nParts = size(data,1);
    t = frames(ii)*ones(nParts,1);
    T = [T; t];
end

