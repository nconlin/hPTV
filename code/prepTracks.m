function tracks = prepTracks(tracks,dt)
% assign naming convention with X as position, T as time, and U as velocity

for ii = 1:numel(tracks)
    if isfield(tracks(ii),'pos')
        tracks(ii).X = tracks(ii).pos'; % make column vector
        tracks(ii).T = tracks(ii).time;
    end
    tracks(ii).U = tracks(ii).U'/dt; % make column vector
end

tracks = rmfield(tracks,'pos');
tracks = rmfield(tracks,'time');