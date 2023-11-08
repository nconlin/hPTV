function [posw,intsw] = newParticleFinder(frm,th,sz,ker)

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

% make sure we have no bad points
good = find(isfinite(xcenters) & isfinite(ycenters));

% fix up the funny coordinate system used by matlab
pos = [ycenters(good), xcenters(good)];

% solve for intensities
ints = NaN(size(pos,1),1);
sigma = NaN(size(pos,1),1);
A = [1 1; 0 1; 1 1];
for ii = 1:size(pos,1)
    
    % format matrix for fitting in x
    b = [z1x(ii); z2(ii); z3x(ii)];
    x = A\b;
    sigma(ii) = (x(1)*-2)^(-1);
    
    % check if signal didnt change (uint8 thing...)
    if abs(sigma(ii)) > 30
        sigma(ii) = 10;
    end

    % solve for x intensity
    intx = exp(x(2))*sqrt(2*pi)*sigma(ii);

    % format matrix for fitting in y
    b = [z1y(ii); z2(ii); z3y(ii)];
    x = A\b;
    sigma(ii) = (x(1)*-2)^(-1);
    
    % check if signal didnt change (uint8 thing...)
    if abs(sigma(ii)) > 30
        sigma(ii) = 10;
    end

    % solve for x intensity
    inty = exp(x(2))*sqrt(2*pi)*sigma(ii);


    % now average x and y
    ints(ii) = 0.5*(intx+ inty);
end

%% average by intensity if they are close together (within sz)
cIdx = clusterdata(pos,'criterion','distance','cutoff',sz);
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

stop = 1;





end



