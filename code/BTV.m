%% CODE TO RUN BUBBBLE TRACKING VELOCIMETRY
% Navigate to a folder where output files can be saved. Run the first section 
% of this script to add the codes to path (may have to copy and paste it to command line for some reason). Then adjust the parameters as
% needed and run the whole thing. 

% clear everything
clear all, close all, clc

%% add codes to path and set directory
addpath('E:\codes\cluster_v2\');
vidDir = 'E:\island_beach\10_9a\';
%% determine number of videos
cd('cam 1\he bubble\')
d = dir;
nVid = length(d)-2;

%%
cd(vidDir)
for ii = 1:4
    fold = ['cam ' num2str(ii) '/bub/'];
    cd(fold);

    d = dir;
    for jj = 3:length(d)
        vidNames{ii,jj-2} = [d(1).folder '\' d(jj).name];
        s = dir(d(jj).name);
        frm(ii,jj-2) = s.bytes;
    end
    cd(vidDir)
end
%% parameters for BTV



% parameter saving
opt = 2;
if exist('params.mat','file')
    disp('params file already exists. would you like to use saved params (1) or set new ones (2)?')
    opt = input('');
end

if opt == 1 % using saved params
    load('params.mat');

elseif opt == 2 % setting new params
    % video files to read
    vidFiles = {'D:\island_beach\10_5\cam 1\he bubble\5M6A1017.mp4', ... 
               'D:\island_beach\10_5\cam 2\he bubble\5M2A0665.mp4', ...
               'D:\island_beach\10_5\cam 3\he bubble\0C5A0503.mp4', ...
               'D:\island_beach\10_5\cam 4\he bubble\0D9A0595.mp4'
                }; % paths to video files
    
    % background params
    window = 500; % window size in frames for each background
    step = 10; % step size in frames to take through video
    
    % center finding params (can be tuned for each video)
    BackgroundType = 'BackgroundMean';
    th = [5 5 5 5]; % intensity threshold
    sz = [5 5 5 5]; % particle size
    
    % 2D tracking params
    maxDisp = 20; % maximum displacement from predicted location (px)
    minLen = 30; % minimum track length in frames
    minDisp = 20; % minimum distance between any two positions in a track 
    
    % track pairing params
    maxDist = 100;
    maxRelVelDiff = 0.4;
    minLenSim = 5;
    
    % offset file
    offsetFile = 'offsets.txt';
    
    % ray conversion
    calibFile = 'E:\island_beach\10_9a\axis_cal_dltCoefs.csv';
    outputFile = 'rays'; % name of file to write outputs (no extension)
    save('params.mat');
else
    error('not a valid option');
end

% compute background (this takes a while)



%% running analysis

runFiles = 11;
times = NaN(length(runFiles),1);
count = 1;
for iter = runFiles

    % go into the folder
    tic
    fold = ['bub_' num2str(iter)];
    cd(fold);
    
    
    % load in the videos
    load('vidFiles.mat')
    load('params_new.mat')

    % backgrounds
    if ~exist('cam4_bgs.mat','file')
        disp('backgrounds')
        window = 1000;
        step = 50;
        bgFiles = {'cam1_bgs.mat','cam2_bgs.mat','cam3_bgs.mat','cam4_bgs.mat'}; % where to save backgrounds
        for ii = 1:4
            V = VideoReader(vidFiles{ii});
            BackgroundComputation_video_window(V,1,V.NumFrames,bgFiles{ii},window,step);
        end
    end

    % identify centers (this is slow) parameters may need to be adjusted
    if ~exist('cam4_centers_new.mat','file')
        centerFiles = {'cam1_centers_new.mat','cam2_centers_new.mat','cam3_centers_new.mat','cam4_centers_new.mat'};
        for ii = 1:4
            disp('centers')
            V = VideoReader(vidFiles{ii});
            CC = comboCenterFinding2D_video_window(V,bgFiles{ii},1:V.NumFrames,th(ii),sz(ii,:),msk_th(ii),centerFiles{ii},BackgroundType);
        end
    end
    
    % synchronization
    if ~exist('offsets.txt','file')
        disp('offsets')
        [offsets_synch,gdness] = synchronizeVideos(vidFiles,1:200,[],cropRect);
        writematrix(offsets_synch,'offsets.txt')
    end

    % save settings
    save('params_new.mat','window','step','BackgroundType','sz','th','offsetFile','calibFile','msk_th','gdness','cropRect')

    % shift for offset 
    offsets = importdata(offsetFile);
    offCenters = {'cam1_centers_off_new.mat','cam2_centers_off_new.mat','cam3_centers_off_new.mat','cam4_centers_off_new.mat'};
    for ii = 1:numel(offCenters)
        load(centerFiles{ii});
        % offsets are negative and tell which frames to discard
        start = abs(offsets(ii))+1;
        CC = CC(start:end);
        nframes = numel(CC);
        save(offCenters{ii},'CC','nframes');
    end
    
    % convert camera points to rays
    [P,V]=Centers2RaysDLT(calibFile,offCenters,outputFile);
    

     % record duration
    cd ..
    iter
    times(count) = toc
    count = count + 1;
end
   
%%
%{
%% convert to tracking format and track in 2D
trkFiles = {'cam1_tracks_2D.mat','cam2_tracks_2D.mat','cam3_tracks_2D.mat','cam4_tracks_2D.mat'};
for ii = 1:numel(centerFiles)
    data = load(centerFiles{ii});
    [X,T] = CCtoTrackData(data.CC);
    [tracks,vtracks] = trackParticles(X,T,maxDisp,minLen,minDisp);
    save(trkFiles{ii},'tracks','vtracks','maxDisp','minLen');
end
%}
%{
%% synchronization test
endFrames = 50:5:200;
offsets_synch = NaN(4,3,length(endFrames));
gdness = NaN(4,3,length(endFrames));
for ii = 1:length(endFrames)
    [offsets_synch(:,:,ii),gdness(:,:,ii)] = synchronizeVideos(vidFiles,1:endFrames(ii),[],cropRect);
end
%}

%% synchronization
[offsets_synch,gdness] = synchronizeVideos(vidFiles,1:200,[],cropRect);
writematrix(offsets_synch,'offsets.txt')
%{
%% pair them up 
for ii = 1:numel(trkFiles)
    load(trkFiles{ii})
    [dualTracks,filterTracks] = pairGlareTracks_v2(vtracks,maxDist,maxRelVelDiff,minLenSim);
    save(trkFiles{ii},'dualTracks','filterTracks','-append')
end
%}
%% convert tracks to centers and save
trackCenters = {'cam1_track_centers_off.mat','cam2_track_centers_off.mat','cam3_track_centers_off.mat','cam4_track_centers_off.mat'};
offsets = importdata(offsetFile);
for ii = 1:numel(trackCenters)
    data = load(trkFiles{ii});
    CC = TrackDatatoCC(data.tracks);
    % offsets are negative and tell which frames to discard
    start = abs(offsets(ii))+1;
    CC = CC(start:end);
    nframes = numel(CC);
    save(trackCenters{ii},'CC','nframes');
end




%% now run STM ... 

% some general sizing considerations
bubDiameter = 0.007; % (m) bubble diameter
xWidth = 10; % (m) total size of domain in x direction
yWidth = 10; 
zWidth = 5;

% approximate number of voxels to use (with a fudge factor)
nX = 1.1*xWidth/bubDiameter
nY = 1.1*yWidth/bubDiameter
nZ = 1.1*zWidth/bubDiameter

% midpoint distance to use
d = 0.8*bubDiameter
%% read in the results

for iter = 8:16

    % go into the folder
    fold = ['bub_' num2str(iter)];
    mkdir(fold)
    cd(fold)

    for pp = 1:4
        vidFiles{pp} = vidNames{pp,iter};
    end
    save('vidFiles.mat','vidFiles')

%% read in the results from matching
matchFile = 'rays_out_cpp.h5';
minFrame = inf;
for ii = 1:4
    load(['cam' num2str(ii) '_centers_off.mat'])
    if numel(CC) < minFrame
        minFrame = numel(CC);
    end
end
frames = 1:minFrame;
[X,T,E] = readSTM(matchFile,frames);

%{
%% merge close together points 

dXtol = 0.01; % distance tolerance to merge together

% find nearest neighbor distances
Xnew = NaN(size(X));
Tnew = NaN(size(T));
Enew = NaN(size(E));

% setup
[frmNums,startInds] = unique(T); % when do frames start and end
endInds = circshift(startInds,-1) - 1;
endInds(end) = size(X,1); 
nFrames = length(frmNums);
fillInds = 1;
startFlag = 1;

for ii = 1:nFrames
    
    % get current positions
    nowInds = startInds(ii):endInds(ii);
    xnow = X(nowInds,:);
    enow = E(nowInds);
    npts = size(xnow,1);
    
    if npts > 1 % do the clustering
    
        % cluster points based on distance
        cIdx = clusterdata(xnow,'criterion','distance','cutoff',dXtol);
        CidxUniq = unique(cIdx);
        xnew = NaN(length(CidxUniq),3);
        enew = NaN(length(CidxUniq),1);
        for jj = 1:length(CidxUniq)
            cnowIdx = find(cIdx == CidxUniq(jj));
            xnew(jj,:) = mean(xnow(cnowIdx,:),1);
            enew(jj) = mean(enow(cnowIdx),1);
        end

    else % just record the point we have

        xnew = xnow;
        enew = enow;

    end

    % record it to the big matrix
    startInd = fillInds(end)+1;
    if startFlag
        startInd = 1;
        startFlag = 0;
    end
    endInd = startInd + size(xnew,1) - 1;
    fillInds = startInd:endInd;
    Xnew(fillInds,:) = xnew;
    Tnew(fillInds) = frmNums(ii);
    Enew(fillInds) = enew;

    
   % report progress
   ii/nFrames
end

% remove extra points
gdInds = find(~isnan(Xnew(:,1)));
Xmerge = Xnew(gdInds,:);
Tmerge = Tnew(gdInds);
Emerge = Enew(gdInds);
%}
%% tracking
% tracking setup 
maxDist = 0.1; % max distance from predicted position (m)
minLen = 5; % min track length in frames
minDisp = 0.05; % min displacement between any two positions (m)

% tracking in 3D (3BE tracker)
[tracks,vtracks] = trackParticles(X,T,maxDist,minLen,minDisp);
save('rawTrackData.mat','tracks','vtracks');


%% yield rate
offCenters = {'cam1_track_centers_off.mat','cam2_track_centers_off.mat','cam3_track_centers_off.mat','cam4_track_centers_off.mat'};
%offCenters = {'cam1_centers_off.mat','cam2_centers_off.mat','cam3_centers_off.mat','cam4_centers_off.mat'};
load('rawTrackData')

avgPartPerFrame = zeros(length(offCenters),1);
avPart = @(S) length(S.X);
for ii = 1:length(offCenters)

    load(offCenters{ii})
    allLens = arrayfun(avPart,CC);
    allLens = allLens(allLens > 0);
    avgPartPerFrame(ii) = mean(allLens);
end

allT = horzcat(tracks.time);
c = unique(allT);
counts = zeros(length(c),1);
for ii = 1:length(c)
    counts(ii) = sum(allT == c(ii));
end

yieldRate = mean(counts)/mean(avgPartPerFrame)
save('yieldrate.mat','yieldRate')

cd ..
end



%% 4d-ptv track
data = [Tmerge Xmerge];
maxDist = 0.1;
minLen = 5;
predict = 1; % flag for predictive tracking
nPrior = 3;  % # of prior of frames to use
conflict = 0; % flag to deal with conflicts
trkFileName = 'exp_test';

[tracks,traj]= track3d_manualfit(' ',trkFileName,data,maxDist,minLen,predict,nPrior,conflict);

%% stitching tracks (only if using 4d-ptv tracker)
dfmax = 5; % maximum number of missing frames tolerated
dxmax = 0.1; % maximum space difference tolerated between expected and real positions
dvmax = 0.3;  %maximum relative velocity difference tolerated
lmin = 15; % minimum length to tru to reconnect a trajectory.
StitchedTraj = stitchTracks(traj,dfmax,[],dxmax,dvmax,lmin);
%}
%% compute yield rate of experiment

%% record video times in text file
for ii = 1:17
    fold = ['bub_' num2str(ii)];
    cd(fold)
    load('vidFiles.mat')
    V = VideoReader(vidFiles{1});
    d = dir(vidFiles{1});
    startTime = datetime(d.date);
    endTime = startTime + seconds(V.Duration*30/120);

    % write it to a file
    fid = fopen('vid_times.txt','wt');
    fprintf(fid,'%23s\n',[startTime endTime]);
    fclose(fid);

    cd ..
end
    


