function [Slong,Slat] = structure_functions(X,U,T,locs,wids)

% compute structure functions given a point cloud over time
% Inputs:
%   X: 3 x nPt array of positions
%   U: 3 x nPt array of velocities
%   T: 1 x nPt array of frame locations
%   locs: 3 x nLocs
%   wids: 3 x nLocs


nLocs = size(locs,2);

for ii = 1:nLocs % loop over locations

    % find positions, velocities, and times in the bin
    hereInds = find( abs(X(1,:) - locs(1,ii)) < wids(1,ii) & ... 
                   abs(X(2,:) - locs(2,ii)) < wids(2,ii) & ...
                   abs(X(3,:) - locs(3,ii)) < wids(3,ii));
    xhere = X(:,hereInds);
    uhere = U(:,hereInds);
    there = T(hereInds);
    
    % get the position and velocity differences
    temp_twoPt = twoPt_diffs(xhere,uhere,there);

    % format for computation and get components
    dX = [vertcat(temp_twoPt.dx)'; vertcat(temp_twoPt.dy)'; vertcat(temp_twoPt.dz)';];
    dU = [vertcat(temp_twoPt.du)'; vertcat(temp_twoPt.dv)'; vertcat(temp_twoPt.dw)';];
    r = sqrt(sum(dX.^2,1));
    dUL = sum(dX.*dU,1)./r; % get component along separation, normalize by distance magnitude


end

Slong = 1;
Slat = 1;

end

function twoPt = twoPt_diffs(X,U,T)
plotit = 0;
% now loop over times
tu = unique(T);
blnk = struct('dx',[],'dy',[],'dz',[],'du',[],'dv',[],'dw',[]);
twoPt = repmat(blnk,[length(tu) 1]);
for jj = 1:length(tu)

    % find which particles are during this time
    nowInds = find(T == tu(jj));

    if numel(nowInds) > 2
        xnow = X(:,nowInds);
        unow = U(:,nowInds);
        nf = size(xnow,2);
    
        % compute their differences
        xf = xnow(1,:); % split into components
        yf = xnow(2,:);
        zf = xnow(3,:); 
    
        % position separations for frame
        dx=reshape(triu(xf-xf',1)',[],1);
        dy=reshape(triu(yf-yf',1)',[],1);
        dz=reshape(triu(zf-zf',1)',[],1);
    
        % now do the same for velocities
        uf = unow(1,:);
        vf = unow(2,:);
        wf = unow(3,:);
    
        % position separations for frame
        du=reshape(triu(uf-uf',1)',[],1);
        dv=reshape(triu(vf-vf',1)',[],1);
        dw=reshape(triu(wf-wf',1)',[],1);
    
        % throw out self points
        bdInds = find(dx == 0); % should be all of them
        dx(bdInds) = [];
        dy(bdInds) = [];
        dz(bdInds) = [];
        du(bdInds) = [];
        dv(bdInds) = [];
        dw(bdInds) = [];
        inds = nchoosek(1:nf,2);

        % store it 
        twoPt(jj).dx = dx;
        twoPt(jj).dy = dy;
        twoPt(jj).dz = dz;
        twoPt(jj).du = du;
        twoPt(jj).dv = dv;
        twoPt(jj).dw = dw;
    end

    % make a figure to check it makes sense
    if plotit
        figure
        quiver3(xf,yf,zf,uf,vf,wf), hold on
        for pp = 1:(nf-1)
            plot3([xf(1) xf(1)+dx(pp)],[yf(1) yf(1)+dy(pp)],[zf(1) zf(1)+dz(pp)],'-k')
        end
        input('keep going?')
    end
    
end

end % function