function [X,T,E] = readSTM(file,frames)
% read data from running c version of stm 

% deal with different input options
if numel(frames) == 1
    frames = 1:frames(1);
elseif numel(frames) == 2
    frames = frames(1):frames(2);
end
nframes = length(frames);

% allocate some space in a buffer
bufSize = 1e5;
Xbuf = NaN(bufSize,3);
Tbuf = NaN(bufSize,1);
Ebuf = NaN(bufSize,1);
bufStart = 1;
X = []; T = []; E = []; % output variables
for ii = 1:nframes

    %% read in data for this frame
    dataPath = ['/frame' num2str(frames(ii)) '_xyze'];
    try
        % grab data and check size
        data = h5read(file,dataPath);
        ndata = size(data,1);
        bufEnd = bufStart + ndata - 1;
    catch % frame not found
        warning(['read failed at frame ' num2str(frames(ii))])

        % put buffer data in the output arrays
        bufEnd = bufEnd - ndata; % last entry
        bufInds = 1:bufEnd;
        X = [X; Xbuf(bufInds,:);];
        T = [T; Tbuf(bufInds);];
        E = [E; Ebuf(bufInds);];
        return; % all done

    end

    %% deal with the data
    if bufEnd > bufSize % buffer size exceeded

        % put the data in the output arrays
        bufEnd = bufEnd - ndata; % last entry
        bufInds = 1:bufEnd;
        X = [X; Xbuf(bufInds,:); data(:,1:3);];
        T = [T; Tbuf(bufInds); frames(ii)*ones(ndata,1)];
        E = [E; Ebuf(bufInds); data(:,4);];

        % reset the buffer
        Xbuf = NaN(bufSize,3);
        Tbuf = NaN(bufSize,1);
        Ebuf = NaN(bufSize,1);
        bufStart = 1;

    elseif ii == nframes % we are done put it directly in the array

        % put the data in the output arrays
        bufEnd = bufEnd - ndata; % last entry
        bufInds = 1:bufEnd;
        X = [X; Xbuf(bufInds,:); data(:,1:3);];
        T = [T; Tbuf(bufInds); frames(ii)*ones(ndata,1)];
        E = [E; Ebuf(bufInds); data(:,4);];
        return; % all done

    else

        % put data in buffer
        bufInds = bufStart:bufEnd;
        Xbuf(bufInds,:) = data(:,1:3);
        Ebuf(bufInds) = data(:,4);
        Tbuf(bufInds) = frames(ii)*ones(ndata,1);

        % update buffer start
        bufStart = bufEnd + 1;
    end

   
end

