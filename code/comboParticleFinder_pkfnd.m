function [posw,intsw] = comboParticleFinder_pkfnd(frm,th,sz1,sz2,ker)
% do the combo particle finder but use pk_fnd instead of cluster data to
% speed up stage one of particle finding
method = 'knnsearch';

% pre-compute logarithims
color_depth = 2^8;
logs = 1:color_depth;
logs = [log(0.0001) log(logs)];

% convert to grayscale if in color
if size(frm,3) ~= 1
    frm = rgb2gray(frm);
end

%% first blur the image (with gaussian of size ker)
frm = imgaussfilt(frm,ker);
s = size(frm);

%% now find local maxima locations (using pk_fnd from 4dptv repo)
% identify the local maxima that are above threshold and within sz(1)  
out = pkfnd(frm,th,sz1);

if ~isempty(out) % we found particles
    x = out(:,2);
    y = out(:,1); % put it back in matlab coordinate system
    
    %% find the positions and intensities
    
    % look up the logarithms of the relevant image intensities
    z1x = logs(frm(sub2ind(s,x-1,y)) + 1)';
    z2 = logs(frm(sub2ind(s,x,y)) + 1)';
    z3x = logs(frm(sub2ind(s,x+1,y)) + 1)';
    
    % compute the centers
    xcenters = -0.5 * (z1x.*(-2*x-1) + z2.*(4*x) + z3x.*(-2*x+1)) ./ ...
        (z1x + z3x - 2*z2);
    
    % do the same for the vertical position
    z1y = logs(frm(sub2ind(s,x,y-1)) + 1)';
    z3y = logs(frm(sub2ind(s,x,y+1)) + 1)';
    ycenters = -0.5 * (z1y.*(-2*y-1) + z2.*(4*y) + z3y.*(-2*y+1)) ./ ...
        (z1y + z3y - 2*z2);
    
    % intensity from averaging nearby points
    e = exp(1);
    ints = mean([z1x z2 z3x z1y z3y]+e,2); % keep log positive definite
    
    % make sure we have no bad points
    good = find(isfinite(xcenters) & isfinite(ycenters) & ints > 0);
    
    % fix up the funny coordinate system used by matlab
    pos = [ycenters(good), xcenters(good)];
    ints = ints(good);
    
    %% average by intensity if they are close together (within sz2)
    if size(pos,1) > 1
        if strcmp(method,'cluster')

            % cluster the data
            cIdx = clusterdata(pos,'criterion','distance','cutoff',sz2);
            CidxUniq = unique(cIdx);
            posw = NaN(length(CidxUniq),2);
            intsw = NaN(length(CidxUniq),1);
            
            
            % loop over each cluster (is there a way around this?)
            for jj = 1:length(CidxUniq)
                
                % find all particles in this cluster
                cnowIdx = find(cIdx == CidxUniq(jj));
                posw(jj,1) = sum(ints(cnowIdx).*pos(cnowIdx,1))/sum(ints(cnowIdx));
                posw(jj,2) = sum(ints(cnowIdx).*pos(cnowIdx,2))/sum(ints(cnowIdx));
                intsw(jj) = mean(ints(cnowIdx));
            
            end

        elseif strcmp(method,'knnsearch')

            knnNum = 20; % # of neighbors to look for
            posw = NaN(size(pos));
            intsw = NaN(size(pos,1),1);
            [idx,dist] = knnsearch(pos,pos,'k',knnNum);
            avail = logical(ones(size(pos,1),1));
            done = 0;
            count = 1;
            while ~done
                if avail(count)
                    inds = idx(count,dist(count,:) < sz2); % find particles within the search radius
                    inds = inds(avail(inds)); % only use if available
                    if length(inds) > 1 % more than one
                        % do weighted average
                        posw(count,1) = sum( ints(inds).*pos(inds,1))/sum(ints(inds));
                        posw(count,2) = sum( ints(inds).*pos(inds,2))/sum(ints(inds));
                        intsw(count) = mean(ints(inds));

                    else % just assign the particle
                        posw(count,:) = pos(count,:);
                        intsw(count) = ints(count);
                    end
                    avail(inds) = 0; % set these particles as unavailable
                end
                count = count + 1;
                done = ~sum(avail);
            end

            % get rid of NaNs
            gdInds = ~isnan(pos(:,1));
            posw = posw(gdInds,:);
            intsw = intsw(gdInds);

        end
    
    else % no particles found
        warning('only one particle')
        posw = pos;
        intsw = ints;
    end % if there are enough particles to cluster
else % if we didnt find particles
    posw = [1 1];
    intsw = 1;

end % function



