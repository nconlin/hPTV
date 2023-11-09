function hndle = plotTracks(tracks,scale)
% plot all tracks in a data set

if nargin < 2
    scale = 0.1;
end
% number of tracks
n = numel(tracks);

% check how data is stored
if ~isfield(tracks,'pos')
    for ii = 1:n
        tracks(ii).pos = [tracks(ii).X];
    end
end

% check dimension
d = min([size(tracks(1).pos,1) size(tracks(1).pos,2)]);

if d == 2
    if ~isfield(tracks,'U')
        for ii = 1:n
            plot(tracks(ii).pos(:,1),tracks(ii).pos(:,2),'o'), hold on, drawnow
        end
    else
        % scatter plot of velocities
        allU = vertcat(tracks.U);
        allUnorm = sqrt( sum(allU.^2,1));
        lowLim = mean(allUnorm) - 2*std(allUnorm);
        highLim = mean(allUnorm) + 2*std(allUnorm);
        for ii = 1:n
            unorm = sqrt( sum( (tracks(ii).U).^2,2 ) );
            scatter(tracks(ii).pos(:,1),tracks(ii).pos(:,2),5,unorm), hold on
        end
        c = colorbar;
        clim([lowLim highLim])
        ylabel(c,'Speed (m/s)')
    end
    xlabel('x (px)'),ylabel('y (px)')
elseif d == 3
    if ~isfield(tracks,'U')
        for ii = 1:n
            plot3(tracks(ii).pos(:,1),tracks(ii).pos(:,2),tracks(ii).pos(:,3),'o'), hold on, drawnow
        end
    else
        
        %{
        % vector field 
        for ii = 1:n
            quiver3(tracks(ii).pos(:,1),tracks(ii).pos(:,2),tracks(ii).pos(:,3),tracks(ii).U(:,1),tracks(ii).U(:,2),tracks(ii).U(:,3),scale), hold on, drawnow
        end
        %}

        % scatter plot of velocities
        allU = horzcat(tracks.U);
        allUnorm = sqrt( sum(allU.^2,1));
        lowLim = mean(allUnorm) - 2*std(allUnorm);
        highLim = mean(allUnorm) + 2*std(allUnorm);
        for ii = 1:n
            unorm = sqrt( sum( (tracks(ii).U).^2,1 ) );
            scatter3(tracks(ii).pos(1,:),tracks(ii).pos(2,:),tracks(ii).pos(3,:),5,unorm), hold on
        end
        c = colorbar;
        clim([lowLim highLim])
        ylabel(c,'Speed (m/s)')


    end
    xlabel('x (m)'),ylabel('y (m)'), zlabel('z (m)')
end
axis equal
hndle = gca;