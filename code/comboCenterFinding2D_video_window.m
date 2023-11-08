function CC = comboCenterFinding2D_video_window(V,BackgroundFile,frames,th,sz,mask_th,centerFile,BackgroundType,channel)
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

%mask_th = 20;
ker = 0.5; % hard coded for now
medFiltSz = [5 5]; % hard coded for now
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
bg_rgb = Background(:,:,:,windowCount);
mask = find(rgb2gray(bg_rgb) > mask_th);
CC = repmat(struct('X',[],'Y',[]),[length(frames) 1]);
%{
try 
    load(centerFile,'kframe','CC','windowCount','frmCount')
    frames = kframe:frames(end);
catch
    save(centerFile,'th','sz','-v7.3','frmCount')
end
%}
for kframe=frames

    if (kframe - windowEdge(windowCount)) >= 0 % we are on the next background window
        windowCount = windowCount + 1;
        try
            bg_rgb = Background(:,:,:,windowCount);
        catch
            bg_rgb = Background(:,:,:,1);
            warning(['Using original background. Could not find requested background # ' num2str(windowCount)])
        end
        mask = find(rgb2gray(bg_rgb) > mask_th);
        try
            save(centerFile,'CC','th','sz','-v7.3','kframe','windowCount','-append')
        catch
            save(centerFile,'CC','th','sz','-v7.3','kframe','windowCount')
        end
    end
    
    Im = read(V,kframe) - bg_rgb;
    if size(Im,3) ~= 1 
        if exist('channel','var')
            Im = Im(:,:,channel);
        else
            Im = rgb2gray(Im);
        end
    end
    
    % mask the image
    Im(mask) = 0;

    % high-pass filter
    Im = Im - medfilt2(Im,medFiltSz);

    % finally find the particles
    pos = comboParticleFinder_pkfnd(Im,th,sz(1),sz(2),ker);
    CC(kframe).X=pos(:,1);
    CC(kframe).Y=pos(:,2);
    frmCount = frmCount + 1;

    % tell how we are doing
    disp(['Frame ' num2str(kframe) ', Particles Found: ' num2str(length(CC(kframe).X))])
end

%% Centers saving into a .mat file
nframes = length(frames);
save(centerFile,"CC",'nframes','th','sz','-v7.3')

end

