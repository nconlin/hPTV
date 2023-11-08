%% uncertainty quantification as outlined in Bhattacharya et al (2020)

%% grab some data

% camera calibration data
calibFile_easyWnd = 'D:\utah_trip_data\9-28-2022\cal_synch_easyWandData.mat';
calibFile_coeffs = 'D:\utah_trip_data\9-28-2022\cal_synch_dltCoefs.csv';
DLTcoeffs = load(calibFile_coeffs);

% camera position data
centerFiles = {'cam1_track_centers_off.mat','cam2_track_centers_off.mat','cam3_track_centers_off.mat','cam4_track_centers_off.mat'};
for ii = 1:length(centerFiles)
    load(centerFiles{ii});
    CCall{ii} = CC;
end

% world position data
stmFile = 'rays_many_out_cpp.h5';
frames = 1:300;
[X,T,E] = readSTM(stmFile,frames);
matchIDs = readSTMmatches(stmFile,frames);

% where to save it
savefile = "uncertainty.mat";

%% estimate the 2D position uncertainty (reprojection)

% disparity vector
d = NaN(size(X,1),2,4);

% loop over the discovered particles
for ii = 1:size(X,1)
    
    % check if we moved to next frame
    Tnow = T(ii);
    if ii == 1
        Tlast = -1;
    end
    if Tnow ~= Tlast
        framCount = 1;
    else
        framCount = framCount + 1;
    end

    % matches for this frame
    matches = matchIDs{Tnow}; 

    % now loop over cameras
    for cam = 1:4

        % project the particle onto the cameras view
        Xc_proj = DLTproj(X(ii,:),DLTcoeffs(:,cam));

        % compare this to the identified position on the camera
        part_id = matches(framCount,cam);
        if ~isnan(part_id) % if part_id = -1 nothing was triangulated
            
            % particle position
            CCcam = CCall{cam};
            Xc_cam(1) = CCcam(Tnow).X(part_id);
            Xc_cam(2) = CCcam(Tnow).Y(part_id);

             % disparity vector
            d(ii,:,cam) = Xc_cam - Xc_proj;

        end
          
    end
    Tlast = Tnow;
end

% averages and variances
d_bar = squeeze(nanmean(d,1));
d_var = squeeze(nanvar(d,1));
Sigma_d = diag(d_var(:));

% save output
if exist(savefile,'file')
    input('woudl you like to overwrite???')
end
save(savefile,'d_bar','d_var','Sigma_d','calibFile_coeffs','stmFile','centerFiles')

%% uncertainty in 2d particle (center finding)

% neglecting this for now ... 
load('params.mat','sz')
temp = NaN(2*4,1);
temp(1:2:end) = sz/(sqrt(2)*3);
temp(2:2:end) = sz/(sqrt(2)*3);
Sigma_det = diag(temp.^2);

%% uncertainty in calibration coefficients
load(calibFile_easyWnd)
close all

% get rmse reprojection error
sigma_d_calx = diag( (easyWandData.dltRMSE).^2 );
sigma_d_caly = diag( (easyWandData.dltRMSE).^2 );

% now the 3d covariance estimate
wand_len = rmoutliers(easyWandData.dDLT);
len_var = var(wand_len)/3;
Sigma_xwCalc = diag([len_var len_var len_var]);
%{
% grab the points
wndpts_wrld = easyWandData.frame5;
badinds = find(isnan(wndpts_wrld(:,1)));
wndpts_wrld(badinds,:) = []; % world coordinates

% camera pts
wndpts_cam = easyWandData.wandPts;
badinds = find(isnan(wndpts_cam(:,1)));
wndpts_cam(badinds,:) = [];
wndpts_cam = [wndpts_cam(:,1:8); wndpts_cam(:,9:16)];

% weed out the outliers
end1 = 1:floor(nCal/2);
end2 = (floor(nCal/2)+1):nCal;
wand_est = sqrt( sum( (wndpts_wrld(end1,:)-wndpts_wrld(end2,:)).^2 , 2) );
len_err = (wand_est - wand_len)./wand_len;

% loop over the usable points
nCal = size(wndpts_cam,1);
use_pts = 1:nCal;
d = NaN(nCal,2,4);
for ii = use_pts

    for jj = 1:4

        % project world point
        Xc_proj = DLTproj(wndpts_wrld(ii,:),DLTcoeffs(:,jj));

        % select cam point
        camind = 2*(jj-1)+1;
        pt_inds = camind:(camind+1);
        Xc_cam = wndpts_cam(ii,pt_inds);
        
        % disparity
        d(ii,:,jj) = Xc_proj - Xc_cam;

    end
end
% averages and variances
d_bar = squeeze(mean(d,1));
d_var = squeeze(var(d,1));
%}

%% now we can back out the covariance in calibration coefficients
ncoeffs = size(DLTcoeffs,1);
Sigma_L = NaN(ncoeffs,ncoeffs,4);

% grab points
wndpts_wrld = easyWandData.frame5;
badinds = find(isnan(wndpts_wrld(:,1)));
wndpts_wrld(badinds,:) = []; % world coordinates
axisPts = size(wndpts_wrld,1)-4 : size(wndpts_wrld,1);
wndpts_wrld( axisPts,:) = [];
nCal = size(wndpts_wrld,1);

% loop over the cameras
for ii = 1:4

    % assemble 2D position covariance matrix
    temp = repmat([sigma_d_calx(ii,ii) sigma_d_caly(ii,ii)],[1 nCal]);
    Sigma_xCal = diag(temp);

    % assemble 3D position covariance matrix and parameter gradient matrix
    A = NaN(2*nCal,ncoeffs);
    Sigma_xwCal = zeros(2*nCal,ncoeffs);
    rowCount = 1;

    for jj = 1:nCal
        
        % this cal point position
        pos = wndpts_wrld(jj,:); 
        
        % put it in the parameter gradient matrix
        [vx,vy] = DLTlin_coeff_vec(pos,DLTcoeffs(:,ii));
        A(rowCount,:) = vx';
        A(rowCount+1,:) = vy';

        % 3D position covariance matrix
        C = DLTlin_pos(pos,DLTcoeffs(:,ii));
        Sigma_xwCal(rowCount,rowCount) = C(1,:)*Sigma_xwCalc*C(1,:)';
        Sigma_xwCal(rowCount+1,rowCount+1) = C(2,:)*Sigma_xwCalc*C(2,:)';

        % increment row count
        rowCount = rowCount + 2;
    end
    
    % now solve the system
    rhs = Sigma_xCal + Sigma_xwCal;
    Sigma_L(:,:,ii) = pinv(A)*rhs*pinv(A)';

end
%% solve for the world coordinate uncertainty

% sum up the reprojection and detection error
Sigma_X = Sigma_d + Sigma_det;

% world coord uncert
Sigma_wrld = NaN(3,3,size(X,1));

% loop over all positions
E = NaN(size(X,1),3); % keep track of error for each position
for ii = 1:size(X,1)

    % grad matrix wrt position
    C = DLTlin_pos(X(ii,:),DLTcoeffs);
    
    % uncertainty from calibration coefficients
    Sigma_a = zeros(8,8);
    rowCount = 1;
    for jj = 1:4
        [vx,vy] = DLTlin_coeff_vec(X(ii,:),DLTcoeffs(:,jj));
        Sigma_a(rowCount,rowCount) = vx'*Sigma_L(:,:,jj)*vx;
        Sigma_a(rowCount+1,rowCount+1) = vy'*Sigma_L(:,:,jj)*vy;
        rowCount = rowCount + 2;
    end
    
    % set up solution and solve with pseudo-inverse`
    rhs = Sigma_a + Sigma_X;
    Sigma_wrld(:,:,ii) = pinv(C)*rhs*pinv(C)';
    E(ii,:) = diag(squeeze(Sigma_wrld(:,:,ii)));   
    
end

save(savefile,'E','Sigma_wrld','Sigma_d','Sigma_det','Sigma_a','-append');

% finally calculate uncertainty 
wrld_std = sqrt(diag(nanmedian(Sigma_wrld,3)))*1e3


figure
lbls = {'$\sigma_x$','$\sigma_y$','$\sigma_z$'};
x_vals = 0:0.01:40;
for ii = 1:3
    pd(ii) = fitdist(squeeze(squeeze(sqrt(Sigma_wrld(ii,ii,:))))*1e3,'Lognormal');
    %histogram(squeeze(squeeze(sqrt(Sigma_wrld(ii,ii,:))))*1e3), hold on
    plot(x_vals,pdf(pd(ii),x_vals),'linewidth',2)
    hold on
end
set(gca,'linewidth',2)
set(gca,'fontsize',20)
legend(lbls,'interpreter','latex')
xlabel('Standard Deviation (mm)','Interpreter','latex')
ylabel('PDF','Interpreter','latex')

%% now we can do velocity uncertainty

sigma_u = sqrt( wrld_std(1)^2 + wrld_std(1)^2 - 0.8*wrld_std(1)*wrld_std(1))*60*1e-3
%% plot the results
figure
hmax = 1;
binEdges = {-20:5:20,-20:5:20};
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
for ii = 1:4
       nexttile
    %subplot(2,2,ii)
    P = hist3(squeeze(d(:,:,ii)),binEdges,'CdataMode','auto','EdgeColor','interp');
    hist3(squeeze(d(:,:,ii)),binEdges,'CdataMode','auto','EdgeColor','interp');
    if max(max(P)) > hmax
        hmax = max(max(P));
    end
    hold on
    %contourf(binEdges{1},binEdges{2},P,'linecolor','none')
    %plot(d_bar(1,ii),d_bar(2,ii),'ok','markerfacecolor','k','markersize',5)
    
    view(2)
    title(['Cam ' num2str(ii)],'interpreter','latex')
    xlim([-20 20])
    ylim([-20 20])
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',12)
 
end
for ii = 1:4
    nexttile(ii)
    caxis([0 0.8*hmax])
end
nexttile(3)
xlabel('$\delta x$ (px)','interpreter','latex','FontSize',15)
nexttile(4)
xlabel('$\delta x$ (px)','interpreter','latex','FontSize',15)
%yticks([])

nexttile(1)
ylabel('$\delta y$ (px)','interpreter','latex','FontSize',15)
%xticks([])
nexttile(3)
ylabel('$\delta y$ (px)','interpreter','latex','FontSize',15)

nexttile(2)
%xticks([])
%yticks([])

hp4 = get(nexttile(4),'Position');
c = colorbar('Position', [hp4(1)+hp4(3)+0.03  hp4(2)+.1  0.03  hp4(2)+hp4(3)*1.5]);
c.TickLabelInterpreter = 'latex';

set(gca,'color','none')



