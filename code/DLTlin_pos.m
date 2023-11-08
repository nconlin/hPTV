function C = DLTlin_pos(pos,all_coeffs)
% get matrix of gradients of DLT function wrt position
% C \in \mathbb{R}^{2 N_{cam} \times 3}

% set up problem
nCam = size(all_coeffs,2);
C = NaN(2*nCam,3);

    for ii = 1:nCam
        % grab this cam calib and pre-compute denominator numerator
        coeffs = all_coeffs(:,ii);
        denom = coeffs(9)*pos(1) + coeffs(10)*pos(2) + coeffs(11)*pos(3) + 1;
        num_x = coeffs(1)*pos(1) + coeffs(2)*pos(2) + coeffs(3)*pos(3) + coeffs(4);
        num_y = coeffs(5)*pos(1) + coeffs(6)*pos(2) + coeffs(7)*pos(3) + coeffs(8);
    
        % indices
        x_ind = 2*(ii-1) + 1;
        y_ind = x_ind + 1;
    
        % x coord
        C(x_ind,1) = (denom*coeffs(1) - num_x*coeffs(9))/denom^2;
        C(x_ind,2) = (denom*coeffs(2) - num_x*coeffs(10))/denom^2;
        C(x_ind,3) = (denom*coeffs(3) - num_x*coeffs(11))/denom^2;
        
        % y coord
        C(y_ind,1) = (denom*coeffs(5) - num_y*coeffs(9))/denom^2;
        C(y_ind,2) = (denom*coeffs(6) - num_y*coeffs(10))/denom^2;
        C(y_ind,3) = (denom*coeffs(7) - num_y*coeffs(11))/denom^2;
    end

end
