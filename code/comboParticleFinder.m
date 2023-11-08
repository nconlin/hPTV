function [posw,intsw] = comboParticleFinder(frm,th,sz1,sz2,ker)

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

%% now find local maxima locations (using code from:
% https://web.stanford.edu/~nto/code/ParticleFinder.m)

% identify the local maxima that are above threshold  
maxes = find(frm >= th & ...
    frm >= circshift(frm,[0 1]) & ...
    frm >= circshift(frm,[0 -1]) & ...
    frm >= circshift(frm,[1 0]) & ...
    frm >= circshift(frm,[-1 0]));

% now turn these into subscripts
s = size(frm);
[x,y] = ind2sub(s, maxes);

% throw out unreliable maxes in the outer ring
good = find(x~=1 & y~=1 & x~=s(1) & y~=s(2));
x = x(good);
y = y(good);

% find the horizontal positions

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

%% merge if they are close together (within sz1)
if size(pos,1) > 1
    cIdx = clusterdata(pos,'criterion','distance','cutoff',sz1);
    CidxUniq = unique(cIdx);
    posmerge = NaN(length(CidxUniq),2);
    intsmerge = NaN(length(CidxUniq),1);
    
    
    % loop over each cluster (is there a way around this?)
    for jj = 1:length(CidxUniq)
        
        % find all particles in this cluster
        cnowIdx = find(cIdx == CidxUniq(jj));
        [~,maxInd] = max(ints(cnowIdx));
        posmerge(jj,:) = pos(cnowIdx(maxInd),:);
        intsmerge(jj) = ints(cnowIdx(maxInd));
    
    end

    %% average by intensity if they are close together (within sz2)
    if size(posmerge,1) > 1
        cIdx = clusterdata(posmerge,'criterion','distance','cutoff',sz2);
        CidxUniq = unique(cIdx);
        posw = NaN(length(CidxUniq),2);
        intsw = NaN(length(CidxUniq),1);
        
        
        % loop over each cluster (is there a way around this?)
        for jj = 1:length(CidxUniq)
            
            % find all particles in this cluster
            cnowIdx = find(cIdx == CidxUniq(jj));
            posw(jj,1) = sum(intsmerge(cnowIdx).*posmerge(cnowIdx,1))/sum(intsmerge(cnowIdx));
            posw(jj,2) = sum(intsmerge(cnowIdx).*posmerge(cnowIdx,2))/sum(intsmerge(cnowIdx));
            intsw(jj) = mean(intsmerge(cnowIdx));
        
        end
    end

else
    warning('only one particle')
    posw = pos;
    intsw = ints;
end

stop = 1;

end



