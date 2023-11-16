function stats = lagrangian_analysis(tracks,h,widths)
% function to compute lagrangian quantities at different heights

% prune tracks
minTrackLen = 30;
minLenSeg = 10;
disp(['throwing out tracks shorter than ' num2str(minTrackLen) ' frames'])
tracks = tracks( vertcat(tracks.len) > minTrackLen);

% put tracks in array
X = horzcat(tracks.X);
U = horzcat(tracks.U);
L = horzcat(tracks.len);

% setup start/stop indices of tracks
startInds = [1 cumsum(L)+1];
startInds(end) = [];
stopInds = startInds - 1;
stopInds(1) = [];
stopInds = [stopInds size(X,2)];
trackInds = NaN(size(X,2),1);
timeInds = NaN(size(X,2),1);
for ii = 1:numel(tracks) % TODO: make this not dumb
    fillInds = startInds(ii):stopInds(ii);
    trackInds(fillInds) = ii*ones(length(fillInds),1); 
    timeInds(fillInds) = 1:length(fillInds);
end

% initialize arrays and struct
blnk = struct('dXbar',[],'dX2',[],'Ubar',[],'dU2',[]);
stats = repmat(blnk,[length(h) 1]);

% loop over heights to find tracks to analyze
for ii = 1:length(h)
    
    % find tracks to compute statistics on 
    allInds = find(abs(X(3,:) - h(ii)) < widths(ii)); % all observations in bin
    trks = trackInds(allInds); % tracks which pass through the bin
    [uniqTracks,uniqTrackInds] = unique(trks,'first'); % each track can only contribute to a bin once

    startInds_h = allInds(uniqTrackInds);
    stopInds_h = stopInds(uniqTracks);
    allLens_h = stopInds_h - startInds_h + 1;
    nTracks_h = length(stopInds_h); % # of tracks for this bin

    % throw out track segments if too short
    disp(['throwing out tracks segments shorter than ' num2str(minLenSeg) ' frames'])
    gdTracks = find(allLens_h > minLenSeg);
    startInds_h = startInds_h(gdTracks);
    stopInds_h = stopInds_h(gdTracks);
    allLens_h = allLens_h(gdTracks);
    maxLen = max(allLens_h);

    % collect position and velocity at this height along all trajectories
    dXarray = NaN(3,maxLen,nTracks_h);
    Uarray = NaN(3,maxLen,nTracks_h);
    for jj = 1:length(startInds_h)
        thisLen = length(startInds_h(jj):stopInds_h(jj));
        dXarray(:,1:thisLen,jj) = X(:,startInds_h(jj):stopInds_h(jj)) - X(:,startInds_h(jj));
        Uarray(:,1:thisLen,jj) = U(:,startInds_h(jj):stopInds_h(jj));
    end

    % average
    dXbar = mean(dXarray,3,'omitnan');
    Ubar = mean(Uarray,3,'omitnan');

    % mean squared displacements
    dX2 = mean( (dXarray - repmat(dXbar,[1 1 size(dXarray,3)])).^2, 3, 'omitnan' ); % < dX^2 >
    dXp = dXarray - repmat(dXbar,[1 1 size(dXarray,3)]);

    % velocity fluctuation correlations
    up = Uarray - repmat(Ubar,[1 1 size(Uarray,3)]);
    R = NaN(3,maxLen);
    Rstd = NaN(3,maxLen);
    S = NaN(3,maxLen);
    F = NaN(3,maxLen);
    sigma = NaN(3,3,maxLen);
    Sstd = NaN(3,maxLen);
    Spdf = cell(3,maxLen);
    Sbins = cell(3,maxLen);
    dXpdf = cell(3,maxLen);
    dXbins = cell(3,maxLen);
    dXstd = NaN(3,maxLen);
    t = 1:maxLen;
    nOrder = 2;
    nSamps = NaN(maxLen,1);
    for jj = 1:maxLen

        % current and time shifted indices
        current = 1:(maxLen-t(jj));
        shift = current + t(jj);

        % cross-correlation
        R(:,jj) = mean( up(:,current,:).*up(:,shift,:),[2 3],'omitnan');
        Rstd(:,jj) = std( up(:,current,:).*up(:,shift,:),0,[2 3],'omitnan');
        nSamps(jj) = mean(sum( ~isnan(up(:,current,:).*up(:,shift,:)),[2 3]),1);


        % structure function
        du = up(:,shift,:) - up(:,current,:);
        dXt = dXp(:,shift,:) - dXp(:,current,:);
        S(:,jj) = mean( du.^nOrder , [2 3], 'omitnan');
       % sigma(:,:,jj) = mean( du*du',[3 2],'omitnan');
        Sstd(:,jj) = std( du.^nOrder ,0, [2 3], 'omitnan');
        F(:,jj) = mean( du.^4,[2 3],'omitnan')./S(:,jj).^2;

        % pdf of velocity increments
        for kk = 1:3
            ducomp = squeeze(squeeze(du(kk,:,:)));
            [Spdf{kk,jj},Sbins{kk,jj}] = histcounts(ducomp(~isnan(ducomp)),'Normalization','pdf');
            Sstd(kk,jj) = std(ducomp(~isnan(ducomp)));
        end

        % pdf of displacements
        for kk = 1:3
            dXcomp = squeeze(squeeze(dXt(kk,:,:)));
            [dXpdf{kk,jj},dXbins{kk,jj}] = histcounts(ducomp(~isnan(dXcomp)),'Normalization','pdf');
            dXstd(kk,jj) = std(ducomp(~isnan(dXcomp)));
        end

    end

    % two particle dispersion

        
    % put it in the struct
    stats(ii).h = h(ii);
    stats(ii).dXbar = dXbar;
    stats(ii).dX2 = dX2;
    stats(ii).Ubar = Ubar;
    stats(ii).R = R;
    stats(ii).Rstd = Rstd;
    stats(ii).S = S;
    stats(ii).Sstd = Sstd;
    stats(ii).Spdf = Spdf;
    stats(ii).Sbins = Sbins;
    stats(ii).nSamps = nSamps;
    stats(ii).F = F;
    stats(ii).uf = var(up,0,[2 3],'omitnan');
    stats(ii).dXpdf = dXpdf;
    stats(ii).dXbins = dXbins;
    stats(ii).dXstd = dXstd;

    % report progress
    disp(['Statistics computed for h = ' num2str(h(ii))])

end
    


