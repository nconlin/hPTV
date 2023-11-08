function [offsets,goodness] = synchronizeVideos(vidFiles,frames,startFrames,cropRect)
% return offsets for synchronizing videos. user selects a region in each
% video to autocorrelate intensities. for now a single offset is used. the
% goodness quantifies how aligned the autocorrelations are over the whole
% video

% initialize intensity array
nv = numel(vidFiles);
intens = cell(nv,1);

% first select regions to correlate (if not assigned)
for ii = 1:length(vidFiles)

    % read in the video
    V = VideoReader(vidFiles{ii});
    f = read(V,1);
    [nx,ny,nz] = size(f);
    
    if ~exist('cropRect','var')
        % ask user for selection
        fig = figure;
        imshow(f)
        [xs,ys] = ginput(4); % just one rectangle for now
        close(fig)
        xl = min(xs); % lower x-bound
        yl = min(ys); % lower y-bound
        xu = max(xs); % upper x-bound
        yu = max(ys); % upper y-bound
    
        % get indices for rectangle
        frmxInds = 1:nx;
        frmyInds = 1:ny;
        autoIndsx = find( frmxInds > xl & frmxInds < xu);
        autoIndsy = find( frmyInds > yl & frmyInds < yu);
    else
        % get indices for rectangle
        autoIndsx = cropRect(ii,1):cropRect(ii,2);
        autoIndsy = cropRect(ii,3):cropRect(ii,4);
    end

    % now compute intensity in all 3 channels
    if ~exist('startFrames','var') || isempty('startFrames')
        startFrames = zeros(nv,1);
    end
    if ~exist('frames','var')
        frames = startFrames + repmat(1:V.NumFrames,[nv 1]);
    end
    tempInt = NaN(length(frames),nz);
    frmCount = 1;
    for jj = 1:length(frames)

        % get the intensity
        frm = read(V,frames(jj));
        for kk = 1:nz
            tempInt(frmCount,kk) = sum(frm(autoIndsy,autoIndsx,kk),'all'); % weird matlab coordinate system (x->y)
        end
        frmCount = frmCount + 1;

    end
    
    % record intensity
    intens{ii} = tempInt;

end


%% now compute autocorrelations

% initialize
figure
subplot(nv-1,1,1)
offsets = NaN(nv,nz);
goodness = NaN(nv,nz);
masterSig_all = intens{1};
masterSig = masterSig_all - mean(masterSig_all,1);

% loop over other cameras
for ii = 2:nv
    
    % grab signal from other camera
    secondSig_all = intens{ii}; 
    secondSig = secondSig_all - mean(secondSig_all,1);

    % autocorrelate
    tempOffset = NaN(nz,1);
    tempGdness = NaN(nz,1);
    for jj = 1:nz
        [r,lags] = xcorr(masterSig(:,jj),secondSig(:,jj),'normalized');
        [tempGdness(jj),lag_ind] = max(r);
        tempOffset(jj) = lags(lag_ind);

        % plot
        subplot(nv-1,1,ii-1)
        plot(lags,r,'-o'), hold on
        
    end
    if ii ~= nv
        xticks([])
    end
    
    % record output
    if all(tempOffset == tempOffset(1))
        offsets(ii,1) = tempOffset(1);
        goodness(ii,1) = mean(tempGdness);
    else
        warning('different channels have different offsets, averaging ...')
        for jj = 1:nz
            offsets(ii,jj) = tempOffset(jj);
            goodness(ii,jj) = tempGdness(jj);
        end
    end
    
    if goodness(ii) < 0.5
        warning(['goodness of synchronization is ' num2str(goodness(ii)) 'for cam ' num2str(ii)])
    end

end

% finally format for return
offsets = round(nanmean(offsets,2));
offsets(1) = 0;
[~,max_ind] = max(offsets);
if max_ind ~= 1 % some cameras triggered before cam 1
    offsets = offsets - offsets(max_ind);
end

end

%% xtra function for background subtraction
function frm_bg_sub = subbg(frm,bg,framnum)

    % first find which window
    windInd = find(framnum < bg.windowEdge,1,'last');
    frm_bg_sub = imsubtract(frm,bg.BackgroundMean(:,:,:,windInd));

end