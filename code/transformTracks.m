function vtracks = transformTracks(vtracks,T,flipz)
% rotate and translate track velocity and positions into streamwise, spanwise, and
% vertical coordinates

%% find the current average velocity vector
vmean = zeros(3,1);
nSamp = 0;
for ii = 1:numel(vtracks)
    v = vtracks(ii).U;
    nSamp = nSamp + size(v,2);
    vmean = vmean + sum(v,2);
end

% get mean direction
vmean = vmean/nSamp;
vmeanDir = vmean/norm(vmean);
if abs(vmeanDir(3)) > .1
    %warning('Mean velocity has large component in vertical ...');
end

% now project to horizontal plane
vmeanPlanar = [vmean(1:2); 0]; % artificially set vertical component of mean to zero
Un = vmeanPlanar/norm(vmeanPlanar);
if abs(Un(3)) > .1
    %warning('Mean velocity has large component in vertical ...');
end


%% find the rotation matrix to orient data: 
if flipz
    oldBasis = [1 0 0;
                0 1 0;
                0 0 1;];
    xnew = Un; znew = -oldBasis(:,3); % flip to be positive up 
    ynew = -cross(xnew,znew); % need to also reflect y for proper rotation
    newBasis = [xnew ynew znew];
    warning('assuming axis calibration: +Z axis points directly down');
else
    oldBasis = [1 0 0;
                0 1 0;
                0 0 1;];
    xnew = Un; znew = oldBasis(:,3); % do not flip to be positive up 
    ynew = cross(xnew,znew); % need to also reflect y for proper rotation
    newBasis = [xnew ynew znew];
    warning('assuming axis calibration: +Z axis points directly up');
end

R = zeros(3,3);
for ii = 1:3
    for jj = 1:3
        R(ii,jj) = newBasis(:,ii)'*oldBasis(:,jj);
    end
end


%% now go through the tracks to rotate them all (position and velocity)
badTracks = zeros([numel(vtracks) 1]);
for ii = 1:numel(vtracks)
    vtracks(ii).U = R*vtracks(ii).U;
    vtracks(ii).X = R*vtracks(ii).X-T;
    zpos = vtracks(ii).X(3,:);
    badInd = find(zpos < 0, 1); % look for one instance of negative z position
    if ~isempty(badInd) % get rid of this track
        badTracks(ii) = 1;
    end
end
vtracks = vtracks(~badTracks);
%{
%% make a picture
figure
quiver3(0,0,0,Un(1),Un(2),Un(3)), hold on
for ii = 1:3
    quiver3(0,0,0,newBasis(1,ii),newBasis(2,ii),newBasis(3,ii))
end
legend('mean','xnew','ynew','znew')
%}
