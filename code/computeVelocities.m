function vtracks = computeVelocities(tracks,filterwidth,fitwidth)
% compute velocities based on convolution with a differentiation kernel
% see Ouellette Thesis

if ~exist('filterwidth','var') 
    filterwidth = 1;
end

if ~exist('fitwidth','var')
    fitwidth = 3*filterwidth; % parameters for the differentiation kernel
end

% check for some special cases
% empirically i find these give best results
if fitwidth == 1
    filterwidth = 0.66;
    warning('correcting filter width');
elseif fitwidth == 2
    filterwidth = 0.9; 
    warning('correcting filter width')
end

% check sizes
ntracks = numel(tracks);
d = size(tracks(1).pos,2);

% find constants
I1 = 0; % this evaluates to zero
I2 = 0.5*filterwidth^2*( sqrt(pi)*filterwidth*erf(fitwidth/filterwidth) - 2*fitwidth*exp(-fitwidth^2/filterwidth^2));
Av = -1/I2;

% define the convolution kernel
vkernel = -fitwidth:fitwidth;
vkernel = Av.*vkernel.*exp(-vkernel.^2./filterwidth^2);

% loop over tracks
vtracks = repmat(struct('len',[],'pos',[],'time',[],'U',[]),ntracks,1);
for ii = 1:ntracks
    U = [];
    for jj = 1:d
        U(:,jj) = conv(tracks(ii).pos(:,jj),vkernel,'valid');
    end
    vtracks(ii).U = U;
    vtracks(ii).len = tracks(ii).len - 2*fitwidth;
    vtracks(ii).pos = tracks(ii).pos(fitwidth+1:end-fitwidth,:); 
    vtracks(ii).time =tracks(ii).time(fitwidth+1:end-fitwidth);
end % for ii=1:numel(tracks)