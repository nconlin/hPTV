function CC = newCenterFinding2D_video_window(V,BackgroundFile,frames,th,sz,centerFile,BackgroundType,channel)
%%% Detect particles position in picture and provides their positions in
%%% px from a video by going frame by frame. Based on 4dptv repo function
%--------------------------------------------------------------------------------
%%% Parameters :
%%%     V                          : video reader object
%%%     bgFile                     : .mat file with backgrounds
%%%     frames                     : frames to analyze
%%%     th                         : threshold
%%%     sz                         : typical size of the particles
%%%     BackgroundType (optional)  : determine which background is substracted to pictures. By defaut is equal to BackgroundMean,
%--------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

ker = 0.5; % hard coded for now
%% setup background
if ~exist('BackgroundType','var')
    BackgroundType="BackgroundMean";
end

try
    load(BackgroundFile,BackgroundType);
    load(BackgroundFile,'windowEdge');
    eval(['Background = ' BackgroundType ';']);
catch
    warning('could not find specified background, choosing mean');
    load(BackgroundFile,'BackgroundMean','windowEdge');
    Background = BackgroundMean;
end

if ~exist('windowEdge','var')
    % background is for every frame
    windowEdge = frames(end)+1;
end

frmCount = 1;
windowCount = 1;
for kframe=frames

    if (kframe - windowEdge(windowCount)) >= 0 % we are on the next background window
        windowCount = windowCount + 1;
    end
    
    Im = read(V,kframe) - Background(:,:,:,windowCount);
    if size(Im,3) ~= 1 
        if exist('channel','var')
            Im = Im(:,:,channel);
        else
            Im = rgb2gray(Im);
        end
    end

    pos = newParticleFinder(Im,th,sz,ker);
    CC(frmCount).X=pos(:,1);
    CC(frmCount).Y=pos(:,2);
    frmCount = frmCount + 1;

    % tell how we are doing
    disp(['Frame ' num2str(kframe) ', Particles Found: ' num2str(length(CC(frmCount-1).X))])
end

%% Centers saving into a .mat file
nframes = length(frames);
save(centerFile,"CC",'nframes','th','sz','-v7.3')

end

