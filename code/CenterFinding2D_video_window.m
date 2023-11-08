function CC = CenterFinding2D_video_window(V,BackgroundFile,frames,th,sz,centerFile,BackgroundType,channel)
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
    Nx = size(Im,2);
    Ny = size(Im,1);
    
    out=pkfnd(Im,th,sz); % Provides intensity maxima positions
    npar = size(out,1);
    
    %% We keep only spots with a gaussian shape
    cnt = 0;
    x = [];
    y = [];
    for j = 1:npar
        Nwidth = 1;
        if (out(j,2)-Nwidth >0)&&(out(j,1)-Nwidth>0)&&(out(j,2)+Nwidth<Ny)&&(out(j,1)+Nwidth<Nx)
            cnt = cnt+1;

            Ip = double(Im(out(j,2)-Nwidth:out(j,2)+Nwidth,out(j,1)-Nwidth:out(j,1)+Nwidth));

            x(end+1) = out(j, 1) + 0.5*log(Ip(2,3)/Ip(2,1))/(log((Ip(2,2)*Ip(2,2))/(Ip(2,1)*Ip(2,3))));
            y(end+1) = out(j, 2) + 0.5*log(Ip(3,2)/Ip(1,2))/(log((Ip(2,2)*Ip(2,2))/(Ip(1,2)*Ip(3,2))));
        end
    end
    CC(frmCount).X=x;
    CC(frmCount).Y=y;
    frmCount = frmCount + 1;

    % tell how we are doing
    disp(['Frame ' num2str(kframe) ', Particles Found: ' num2str(length(x))])
end

%% Centers saving into a .mat file
nframes = length(frames);
save(centerFile,"CC",'nframes','th','sz','-v7.3')

end

