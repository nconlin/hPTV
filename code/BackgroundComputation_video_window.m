function BackgroundComputation_video_window(V,StartFrame,EndFrame,saveFile,window,Step)
%%% Compute the background for the pictures of the camera NumCam in the
%%% experiment ManipName
%%% Between picture StartFrame and EndFrame ( every Step frame) it takes
%%% the maximal/minimal/mean intensity for each pixel.
%----------------------------------------------------------------------------
%%% Parameters : 
%%%     V            : video reader object           
%%%     StartFrame   : number of first frame
%%%     EndFrame     : number of the last frame
%%%     saveFile     : path to save the results
%%%     window       : width of window for moving background
%%%     Step         : step between 2 frame taken for computation
%------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% If Step is not define take a step such as the average will be done over
% 1000 pictures.
if (~ exist("Step","var") && (EndFrame-StartFrame)>=1000)
    Step = floor((EndFrame-StartFrame)/1000);
elseif (~ exist("Step","var") && (EndFrame-StartFrame)<1000)
    Step = 1;
end

% allocate space for backgrounds
I = read(V,StartFrame);
dims = size(I);
nBGs = floor((EndFrame-StartFrame+1)/window);
BackgroundMin = zeros([dims(1:2) nBGs]);
BackgroundMax = zeros([dims(1:2) nBGs]);
BackgroundMean = zeros([dims nBGs]); % min and max are non-colored

totsamps = EndFrame - StartFrame + 1;
samples = 1;
windowCount = 1;
windowEdge = StartFrame + window:window:EndFrame;
windowEdge = [windowEdge EndFrame + 1]; % append this so the last window can be checked


% main loop
h = waitbar(0,'Calculating Background ...');
for kframe=StartFrame:Step:EndFrame

    if kframe == StartFrame % just starting out

        I = read(V,kframe);
        BackgroundMintemp = rgb2gray(I);
        BackgroundMaxtemp = rgb2gray(I);
        BackgroundMeantemp = uint32(I);
        disp('window 1')
        disp(num2str(kframe))
        

    elseif (kframe - windowEdge(windowCount)) >= 0 % we stepped to the next window -> start a new one

        % calculate the mean for previous window and record the others
        BackgroundMin(:,:,windowCount) = BackgroundMintemp;
        BackgroundMax(:,:,windowCount) = BackgroundMaxtemp;
        BackgroundMean(:,:,:,windowCount) = BackgroundMeantemp/samples;
        
        windowCount = windowCount + 1; % increment the window
        samples = 1; % reset the sample count 

        % now start reading in the new window
        I = read(V,kframe);
        BackgroundMintemp = rgb2gray(I);
        BackgroundMaxtemp = rgb2gray(I);
        BackgroundMeantemp = uint32(I);

        disp(['Window ' num2str(windowCount)]);
        disp(num2str(kframe));
        
        
    else % nothing special about this frame
        
        I = read(V,kframe);
        BackgroundMaxtemp = max(cat(3,BackgroundMaxtemp,rgb2gray(I)),[],3);
        BackgroundMeantemp = BackgroundMeantemp + uint32(I);
        BackgroundMintemp = min(cat(3,BackgroundMintemp,rgb2gray(I)),[],3);
        samples = samples + 1;
        disp(num2str(kframe))

    end    
    waitbar(kframe/totsamps,h)
end

% assign the final bg 
BackgroundMin(:,:,windowCount) = BackgroundMintemp;
BackgroundMax(:,:,windowCount) = BackgroundMaxtemp;
BackgroundMean(:,:,:,windowCount) = BackgroundMeantemp/samples;

close(h);
BackgroundMean = uint8(BackgroundMean);
save(saveFile,'BackgroundMax','BackgroundMean','BackgroundMin','windowEdge','-v7.3')


