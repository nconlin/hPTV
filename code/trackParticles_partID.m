function [tracks,vtracks] = trackParticles_partID(X,T,maxDisp,minLen,minDisp,matchIDs)
% track particles in arbitrary dimensions with 3 frame best estimate and
% optional 3 frame initialization. Based on PredictiveTracker of Oullete
% Lab

% Inputs:
% X particle positions 
% T time vector of particles
%       eg X(T(1),:) = position of particle 1 and time 1
% maxDisp maximum distance from predicted position
% minLen minimum track length
% minDisp minimum total displacement
% matchIDs matching IDs from particle traingulation
%   cell array where each cell is an array of ids. the row is the particle
%   number (in 3D) and the each column is the particle identifier for each
%   camera

if nargin < 6
    yesPlot = 0;
end
if nargin < 4
    minLen = 5;
end

predictFirstFrame = 0; % TODO make this work
kinPredict = 1; % flag for kinematic prediction

% check for bad input data
badInds = isnan(X(:,1)) | isnan(X(:,2)) | isnan(T);
if ~isempty(badInds)
    warning('Nans found in input data... pruning')
end
X(badInds,:) = [];
T(badInds) = [];
matchIDs(badInds,:) = [];

% sort out times and dimension of data
[frmNums,startInds] = unique(T); % when do frames start and end
endInds = circshift(startInds,-1) - 1;
endInds(end) = size(X,1); 
nFrames = length(frmNums);
d = size(X,2); % dimension of data



% initialize particle tracks
blankTrack = struct('len',[],'pos',[],'time',[],'part_ids',[]);

for k = 1:nFrames
    
    % current frame indices
    nowInds = startInds(k):endInds(k);
    Np = length(nowInds); % particles in this frame

    if Np > 0

        x = X(nowInds,:); % get this frames positions
        partIDs = matchIDs(nowInds,:);
        
        if k == 1 && predictFirstFrame % use 3 frame initialization

            % grab position for next 2 frames
            pos0 = x;
            nextInds = startInds(k+1):endInds(k+1);
            pos1 = X(nextInds,:); 
            nextNextInds = startInds(k+2):endInds(k+2);
            pos2 = X(nextNextInds,:); 

            links = three_frame_init(pos0,pos1,pos2,maxDisp,1); % get the links
            
            % --- now we can start tracks -------
            nLinks = sum(links ~= 0);
            trkCount = 1;
            tracks = repmat(blankTrack,[nLinks 1]);
            for ii = 1:Np
                if links(ii) ~= 0 
                    % this track found a match
                    tracks(trkCount).pos(1,:) = pos0(ii,:);
                    tracks(trkCount).pos(2,:) = pos1(links(ii),:);
                    tracks(trkCount).len = 2;
                    tracks(trkCount).time = [frmNums(1) frmNums(2)];
                    trkCount = trkCount + 1;
                end % if links(ii) ~= 0
            end % for loop particles in initial frame

            % all these tracks are active for now
            active = 1:nLinks;
            Nactive = nLinks;
            disp('tracks started with 3 frame init')

            continue; % done with this frame
           


        elseif k == 2 && predictFirstFrame 
             disp('skipping frame 2')
            continue; % already have tracks so skip


        elseif k == 1 && ~predictFirstFrame

            % just initialize every track
            tracks = repmat(blankTrack,[Np 1]);
            active = 1:Np; % every track is active ( for next frame )
            Nactive = length(active);
            
            % init all tracks
            for ii = 1:Np
                tracks(ii).len = 1;
                tracks(ii).pos = x(ii,:);
                tracks(ii).time = frmNums(k);
                tracks(ii).part_ids = partIDs(ii,:);
            end
           
            % report it out
            disp(['Frame ' num2str(frmNums(k)) ' of ' num2str(nFrames) ': ' ...
            num2str(Np,'%.0f') ' particles, ' ...
            num2str(Np,'%.0f') ' active, ' ...
            num2str(Np,'%.0f') ' new, ' ...
            num2str(1,'%.0f') ' unmatched, ' ...
            num2str(Np,'%.0f') ' total tracks.'])

            % okay we are done with this frame
            continue;

        end % if for first frame predictions/initialization
        
        
        % ---------- do a kinematic prediction of active tracks!! -------------
        pos_prev = NaN(Nactive,d); % positions of two previous frames
        pos_2prev = NaN(Nactive,d);
        
        for ii = 1:length(active)
            pos_prev(ii,:) = tracks(active(ii)).pos(end,:);
            if tracks(active(ii)).len > 1 && kinPredict
                pos_2prev(ii,:) = tracks(active(ii)).pos(end-1,:);
            else
                pos_2prev(ii,:) = tracks(active(ii)).pos(end,:);
            end
        end
        positionGuess = 2*pos_prev - pos_2prev;
       
        
        
        %{
        % grab position data from active tracks
        activeLengths = vertcat(tracks(active).len);
        positionGuess = NaN(Nactive,d);

        % ------- predictions with model based on track length ----------
        for ii = 1:length(active)
            if activeLengths(ii) == 1 % do nearest neighbor ( x_n = x_{n-1})
                positionGuess(ii,:) = tracks(active(ii)).pos(end,:);
            elseif activeLengths(ii) >= 2 % extrapolate from velocity ( x_n = x_{n-1} + x_{n-1} - x_{n-2} )
                positionGuess(ii,:) = 2*tracks(active(ii)).pos(end,:) - tracks(active(ii)).pos(end-1,:);
            %elseif activeLengths(ii) > 2 % extrapolate from velocity and acceleration ( x_n = x_{n-1} + x_{n-1} - x_{n-2} + 0.5*(x_{n-1} - 2*x_{n-2} + x_{n-3})
            %    positionGuess(ii,:) = 2.5*tracks(active(ii)).pos(end,:) - 2*tracks(active(ii)).pos(end-1,:) + 0.5*tracks(active(ii)).pos(end-2,:);
            end
        end
        %}
        
        % ----- evaluate kinematic predictions ---------
        minCost = NaN(Nactive,1); % min cost for each track
        links = zeros(Nactive,1); % init with no particles linked
        dist_est2pos = NaN(Nactive,Np); % distance from guess(i) to particle(j)
        for ii = 1:Nactive
    
            % compute costs
            dist_est2pos(ii,:) = (sum( (positionGuess(ii,:) - x).^2 , 2))'; 
            minCost(ii) = min(dist_est2pos(ii,:));
            
            if minCost(ii) > maxDisp^2 % min cost still too big
                continue; % done with this track, couldnt find a match
            end
    
            % find what the best match is
            bestMatch = find(minCost(ii) == dist_est2pos(ii,:));
            if length(bestMatch) > 1 % uh oh we have a tie
                continue; % just give up for now (Ouellette et al. 2006)
            end
    
            % try to link this track to a new particle
            otherMatchInd = find( links == bestMatch, 1); % find if one other track is already matched to this one
            if ~isempty(otherMatchInd) % this particle has been linked already
                if minCost(otherMatchInd) > minCost(ii) % check if this new link is lower cost
                    links(otherMatchInd) = 0; % then unlink other track
                else
                    links(ii) = 0; % cant link this particle
                    continue; % jump out of loop over track linking
                end
            end
            try
                links(ii) = bestMatch;
            catch
                disp('encountered an error: likely bad input data ... ')
            end
        end % loop over active tracks
    
        % --- now we can append new matches to current tracks -------
        matched = zeros(Np,1);
        for ii = 1:Nactive
            if links(ii) ~= 0 
                % this track found a match
                tracks(active(ii)).pos(end+1,:) = x(links(ii),:);
                tracks(active(ii)).len = tracks(active(ii)).len + 1;
                tracks(active(ii)).time(end+1) = frmNums(k);
                tracks(active(ii)).part_ids(end+1,:) = partIDs(links(ii),:);
                
                % note down that it was matched
                matched(links(ii)) = 1;
            end % if links(ii) ~= 0
        end % for loop over active tracks
    
        % retain which tracks are still active
        active = active(links~=0);
    
        % --- now start new tracks for those unmatched ---- 
        unmatched = find(matched == 0);
        newtracks = repmat(blankTrack,[numel(unmatched) 1]);
        for ii = 1:numel(unmatched)
            newtracks(ii).len=1;
            newtracks(ii).pos = x(unmatched(ii),:);
            newtracks(ii).time = frmNums(k);
            newtracks(ii).part_ids = partIDs(unmatched(ii),:);
        end
        
    else % if there were no particles this frame
        active = [];
        newtracks = [];
        unmatched = [];
        Nactive = 0;
    end % if wrapping for no particles

    % now append progress if tracks already exist
    if exist('tracks','var')
        active = [active (numel(tracks)+1): ...
            (numel(tracks)+numel(newtracks))];
        if ~isempty(newtracks)
            tracks = [tracks ; newtracks];
        end
        Nactive = numel(active);
        links=[];
    end

    disp(['Frame ' num2str(frmNums(k)) ' of ' num2str(nFrames) ': ' ...
        num2str(Np,'%.0f') ' particles, ' ...
        num2str(Nactive,'%.0f') ' active, ' ...
        num2str(numel(unmatched),'%.0f') ' new, ' ...
        num2str(sum(links==0),'%.0f') ' unmatched, ' ...
        num2str(numel(tracks),'%.0f') ' total tracks.'])

end % for loop over frames


% -- pruning ---
disp(['pruning for tracks less than ' num2str(minLen)])
good = vertcat(tracks.len) > minLen;
tracks = tracks(good);
disp([ 'now we have ' num2str(numel(tracks)) ' tracks'])
disp(['pruning for tracks with less than ' num2str(minDisp) ' total displacement'])
totDisp = NaN([numel(tracks) 1]);
for ii = 1:numel(tracks)
    totDisp(ii) = sum( sqrt( sum( (tracks(ii).pos(2:end,:)-tracks(ii).pos(1:end-1,:)).^2,2) ) );
end
good = totDisp > minDisp;
tracks = tracks(good);
disp([ 'now we have ' num2str(numel(tracks)) ' tracks'])


end % function

