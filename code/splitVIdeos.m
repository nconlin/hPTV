%% split up videos script

load('vidNames.mat')
splitInds = 4:7;
splitLen = 10000;

for ii = 1:length(splitInds)

    for jj = 1:4
        
        % load in the video
        curVid = vidNames{jj,splitInds(ii)};
        V = VideoReader(curVid);
        nf = V.NumFrames;

        % setup output
        [fold,name,ext] = fileparts(curVid);



    end
end