function [P,V]=findRaysDLT(calib,x_px,y_px,f,Npoints)
%% calib : calibration data for this camera (DLT coefficients)
%% x_px  : x coordinates in px,
%% y_px  : y coordinates in px,
%% f     : optional focal length input to build interpolation points   

if nargin < 5
    Npoints = 3;
end

% number of particles
Npart = numel(x_px);
Zpts = f*(1:Npoints);

% 3D points from inverting calibration
XYZ = zeros(Npoints,3,Npart);

for ii = 1:Npoints

    % invert the DLT calibration at each point along Z
    for jj = 1:Npart
        if ~isnan(x_px(jj)) && ~isnan(y_px(jj))
            [A,b] = getAmatrix_and_bvec(calib,x_px(jj),y_px(jj),Zpts(ii));
            sol = A\b;
            XYZ(ii,1,jj)=sol(1);
            XYZ(ii,2,jj)=sol(2);
            XYZ(ii,3,jj)=Zpts(ii);
        else
            XYZ(ii,1,jj) = NaN;
            XYZ(ii,2,jj) = NaN;
            XYZ(ii,3,jj) = NaN;
        end
    end
end
[P, V] = fit3Dline(XYZ);
end
%% function for formatting DLT system
function [A,b] = getAmatrix_and_bvec(dlt,x_px,y_px,z_wrld)
% format DLT inversion into linear system

A = [dlt(1)-x_px*dlt(9) dlt(2)-x_px*dlt(10);
     dlt(5)-y_px*dlt(9) dlt(6)-y_px*dlt(10)];

b = [z_wrld*(x_px*dlt(11) - dlt(3)) + x_px - dlt(4);
     z_wrld*(y_px*dlt(11) - dlt(7)) + y_px - dlt(8)];

end