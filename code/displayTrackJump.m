%% display track jumps in 2D

% set up
trkFile = 'cam1_tracks_2D.mat';
centerFile = 'cam1_centers.mat';
load('vidFiles')
vidFile = vidFiles{1};
bgFile = 'cam1_bgs.mat';
load(bgFile)
viewWid = 50;
viewHt = 20;


%% load in data
load(trkFile)
load(centerFile)
V = VideoReader(vidFile);

%% choose a track
trkPick = 2;
[~,idx] = sort(vertcat(tracks.len),'descend');

%% find a jump
tjump = 255;
bgInd = find(tjump<windowEdge);
bgInd = bgInd(end);

%% display images
H = VideoWriter('example_track_jump.avi');
H.Quality = 95;
open(H)

twid = -35:35;
tstart = vtracks(idx(trkPick)).time(1);
posRef = [vtracks(idx(trkPick)).pos(tjump-tstart,1) vtracks(idx(trkPick)).pos(tjump-tstart,2)];
done = 0;
now = tjump;
inc = 0;

th = 5;
sz = 12;
ker = 0.5;
%%
while ~done

    % read in video frame
    now = now + inc;
    Im = rgb2gray(imsubtract(read(V,now),BackgroundMean(:,:,:,bgInd)));
    imshow(Im), hold on

    % all detected particles
    %plot(CC(now).X,CC(now).Y,'x')

    % now the track
   % plot(vtracks(idx(trkPick)).pos(:,1),vtracks(idx(trkPick)).pos(:,2),'-o')

    % zoom in
    xlim([posRef(1)-viewWid posRef(1)+viewWid])
    ylim([posRef(2)-viewHt posRef(2)+viewHt])

    % old particle finder
    Nx = size(Im,2);
    Ny = size(Im,1);
    out=pkfnd(Im,th,sz); % Provides intensity maxima positions
    npar = size(out,1);
    
    % keep only spots with a gaussian shape
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
    b4 = plot(x,y,'o'); hold on

    % use the new particle finder!!!
    sz1 = 3;
    sz2 = 10;
    [newPos,intsw] = comboParticleFinder(Im,th,sz1,sz2,ker);
    aft = plot(newPos(:,1),newPos(:,2),'x');
    legend([b4 aft],'Before','After')


    % wait for input
    drawnow
    h = input('go back (1), go forward (2), done (3)');
    hold off
    if h == 1
        inc = -1;
    elseif h == 2
        inc = 1;
    else
        done = 1;
    end
    vidframe = getframe(gcf);
    writeVideo(H,vidframe);
end
%% 
close(H)
