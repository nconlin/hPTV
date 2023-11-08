function matches = readSTMmatches(file,frames)
% read data from running c version of stm 
nCam = 4;

if numel(frames) == 1
    frames = 1:frames(1);
elseif numel(frames) == 2
    frames = frames(1):frames(2);
end
nframes = length(frames);
matches = cell(nframes,1);

for ii = 1:nframes

    % read in data for this frame
    dataPath = ['/frame' num2str(frames(ii)) '_camrayids'];
    data = h5read(file,dataPath);
    
    % assign cam1-4 ray ids
    matchnow = NaN(size(data,1),nCam);
    for jj = 1:size(data,1) % loop over particles in this frame

        % grab data and prune
        if size(data,2) == 8 % four cameras
            cams = data(jj,[1 3 5 7]);
            rayids = data(jj,[2 4 6 8]);
            badinds = find(cams < 1);
            cams(badinds) = [];
            rayids(badinds) = [];
        elseif size(data,2) ==  6 % three cameras
            cams = data(jj,[1 3 5]);
            rayids = data(jj,[2 4 6]);
            badinds = find(cams < 1);
            cams(badinds) = [];
            rayids(badinds) = [];
            %cams = [cams NaN];
            %rayids = [rayids NaN];
        end
        
        % assign output
        matchnow(jj,cams) = rayids;

    end
  
    matches{ii} = matchnow;

end

