function [P,V]=Centers2RaysDLT(calibFile,centerFiles,outputFile,f)
%   For a given DLT calibration and center files with
%   pixel positions, determines the corresponding rays of light.
%   camID is an array containing the IDs of the cameras center files to be
%   processed (for instance camID=[0 1 2] for a 3 cameras tracking
%   experiment). f is the focal length in world units used as a scale to build points
%   for ray fitting
%----------------------------------------------------------------------------------------
% INPUT Parameters:
%   CalibFile   : Path of the calibration file
%   centerFiles : cell array of paths to particle centers
%   f           : optional input, focal length in world units
% ------------------------------------------------------------------------------------------

% If focal length does not exist set a default value
if ~exist('f','var')
    f = 1; % set to 1 in world units
end
vis = 0;
% set up camera IDs
camID = 1:numel(centerFiles);
%% Folders, Files and Loading
% Load calibration file (contains DLT coefficients along the columns)
fprintf("Loading calibration file... \n")
try
    calib = load(calibFile);
catch
    calib = importdata(calibFile);
    warning('Couldnt load at mat file... Loading as csv')
end

if size(calib) ~= [11 numel(camID)]
    error('calibration should be 11 x nCam')
end

% Open results ray file
fileRays=[outputFile '.dat'];
fid=fopen(fileRays,'w');

% Loop over cameras
disp(['Starting loop over ' num2str(camID(end)) ' cameras...'])
for kcam=1:numel(camID)

    % grab this cameras calibration
    calibNcam=calib(:,kcam);
    
    % load centers for camera camID(kcam)
    fileCenters=centerFiles{kcam};
    [CC,nframes] = readCentersMAT(fileCenters); 
    if kcam == 1
        minFrames = nframes;
    elseif nframes < minFrames
        minFrames = nframes;
    end
    % loop over frames
    for k=1:nframes
        if rem(k,100)==0
            fprintf("cam %d frame %d...\n",kcam,k)
        end
        % convert pixel coordinates into rays of light using the
        % calibration
        
        if isempty(CC(k).X) || isempty(CC(k).Y) || all(isnan(CC(k).Y)) || all(isnan(CC(k).X)) % no particles from this camera on this frame
            fprintf("no particles in camera %d frame %d...\n",kcam,k)
            data(k).P = [1 1 1];
            data(k).V = [1 1 1];
            data(k).rayID = 1;
        else
            [P,V]=findRaysDLT(calibNcam,CC(k).X',CC(k).Y',f);
            rayID=find((~isnan(P(:,1)))&(~isempty(P(:,1))));
            data(k).P=P(rayID,:);
            data(k).V=V(rayID,:);
            data(k).rayID=rayID;
        end
    end
    datacam(kcam).data=data;
end
fprintf("Writing data in matlab file in progress...")
% Writing data in matlab file
save([outputFile '.mat'],'datacam','-v7.3')

%% write results in file

fprintf("Writing data in process...\n")
for kframe=1:minFrames
    Nrays=0;
    for kcam=1:numel(camID)
        Nrays = Nrays + numel(datacam(kcam).data(kframe).rayID);
    end
    fwrite(fid,Nrays,'uint32');
    
    for kcam=1:numel(camID)
        for kray=1:numel(datacam(kcam).data(kframe).rayID)
            fwrite(fid,camID(kcam),'uint8');
            fwrite(fid,datacam(kcam).data(kframe).rayID(kray),'uint16');
            fwrite(fid,datacam(kcam).data(kframe).P(kray,:),'float32');
            fwrite(fid,datacam(kcam).data(kframe).V(kray,:),'float32');
        end
    end
end

fclose(fid);

%% visualize if requested
if vis
    % grab some rays
    figure
    for ii = 1:nframes
        check = 0;
        for jj = 1:numel(datacam)
            datnow = datacam(jj).data; % select camera
            datnow = datnow(ii); % select frame
            if ~isempty(datnow.P) % we have particles
                check = 1;
                for kk = 1:size(datnow.P,1)
                    plotRay3D(datnow.P(kk,:),datnow.V(kk,:)), hold on
                end
                
            end
        end
        if check
            title(num2str(ii))
            input('press enter to continue'), hold off
        end
        
    end
end
