function [E,wrld_std] = uncert_quant_tracks(tracks,CCall,DLTcoeffs)
% function implementation of uncertainty quantification as outlined in Bhattacharya et al (2020)
% computes uncertainty at all tracked points
plotResults = 1;

% world position data from tracks
X = vertcat(tracks.pos); % positions
T = horzcat(tracks.time); % times
matches = vertcat(tracks.part_ids); % particle ids

%% estimate the 2D position uncertainty (reprojection)

% disparity vector
d = NaN(size(X,1),2,4);

% loop over the discovered particles
for ii = 1:size(X,1)
    
    % now loop over cameras
    for cam = 1:4

        % project the particle onto the cameras view
        Xc_proj = DLTproj(X(ii,:),DLTcoeffs(:,cam));

        % compare this to the identified position on the camera
        part_id = matches(ii,cam);
        if ~isnan(part_id) % if part_id = -1 nothing was triangulated
            
            % particle position
            CCcam = CCall{cam};
            Xc_cam(1) = CCcam(T(ii)).X(part_id);
            Xc_cam(2) = CCcam(T(ii)).Y(part_id);

             % disparity vector
            d(ii,:,cam) = Xc_cam - Xc_proj;

        end
          
    end
end

% averages and variances
d_bar = squeeze(nanmean(d,1));
d_var = squeeze(nanvar(d,1));
Sigma_d = diag(d_var(:));


%% solve for the world coordinate uncertainty

% sum up the reprojection and detection error
Sigma_X = Sigma_d;% + Sigma_det;

% world coord uncert
Sigma_wrld = NaN(3,3,size(X,1));

% loop over all positions
E = NaN(size(X,1),3); % keep track of error for each position
for ii = 1:size(X,1)

    % grad matrix wrt position
    C = DLTlin_pos(X(ii,:),DLTcoeffs);
    
    %{
    % uncertainty from calibration coefficients
    Sigma_a = zeros(8,8);
    rowCount = 1;
    for jj = 1:4
        [vx,vy] = DLTlin_coeff_vec(X(ii,:),DLTcoeffs(:,jj));
        Sigma_a(rowCount,rowCount) = vx'*Sigma_L(:,:,jj)*vx;
        Sigma_a(rowCount+1,rowCount+1) = vy'*Sigma_L(:,:,jj)*vy;
        rowCount = rowCount + 2;
    end
    %}
    % set up solution and solve with pseudo-inverse`
    %rhs = Sigma_a + Sigma_X;
    rhs = Sigma_X;
    Sigma_wrld(:,:,ii) = pinv(C)*rhs*pinv(C)';
    E(ii,:) = diag(squeeze(Sigma_wrld(:,:,ii))); % variances in each direction
    
end

%% finally calculate uncertainty 
wrld_std = sqrt(diag(nanmedian(Sigma_wrld,3)))

if plotResults
    %% plot the results
    f = figure;
    f.Units = 'inches';
    f.Position = [0 0 5 5];
    hmax = 1;
    edgeVal = 4;
    binEdges = {-edgeVal:0.5:edgeVal,-edgeVal:0.5:edgeVal};
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
        xlim([-edgeVal edgeVal])
        ylim([-edgeVal edgeVal])
        set(gca,'linewidth',1.5)
        set(gca,'fontsize',12)
     
    end
    for ii = 1:4
        nexttile(ii)
        caxis([0 0.8*hmax])
    end
    nexttile(3)
    xlabel('$\delta x_c$ (px)','interpreter','latex','FontSize',15)
    nexttile(4)
    xlabel('$\delta x_c$ (px)','interpreter','latex','FontSize',15)
    %yticks([])
    
    nexttile(1)
    ylabel('$\delta y_c$ (px)','interpreter','latex','FontSize',15)
    %xticks([])
    nexttile(3)
    ylabel('$\delta y_c$ (px)','interpreter','latex','FontSize',15)
    
    nexttile(2)
    %xticklabels([])
    %yticks([])
    
    hp4 = get(nexttile(4),'Position');
    c = colorbar('Position', [hp4(1)+hp4(3)+0.03  hp4(2)+.1  0.03  hp4(2)+hp4(3)*1.5]);
    c.TickLabelInterpreter = 'latex';
    
    set(gca,'color','none')

    %% 3D position uncertainty
    figure
    lbls = {'$\sigma_x$','$\sigma_y$','$\sigma_z$'};
    x_vals = 0:0.01:40;
    for ii = 1:3
        % fit a pdf
       % pd(ii) = fitdist(squeeze(squeeze(sqrt(Sigma_wrld(ii,ii,:))))*1e3,'Lognormal');
        %plot(x_vals,pdf(pd(ii),x_vals),'linewidth',2)
    
        % true histogram
        [cnts,edges] = histcounts(squeeze(squeeze(sqrt(Sigma_wrld(ii,ii,:))))*1e3,50,'normalization','pdf');
        plot(edges(2:end),cnts,'-o','linewidth',2,'markersize',2)
        
        hold on
    end
    set(gca,'linewidth',2)
    set(gca,'fontsize',20)
    legend(lbls,'interpreter','latex')
    xlabel('Standard Deviation (mm)','Interpreter','latex')
    ylabel('PDF','Interpreter','latex')

end
