function C = DLTlin_pos_num(pos,all_coeffs)
% get matrix of gradients of DLT function wrt position
% C \in \mathbb{R}^{2 N_{cam} \times 3}
% calculate numerically to validate 


% set up problem
nCam = size(all_coeffs,2);
C = NaN(2*nCam,3);
dx = [0.01 0 0];
dy = [0 0.01 0];
dz = [0 0 0.01];

    for ii = 1:nCam

        % grab cam cal
        coeffs = all_coeffs(:,ii);
    
        % perturb in each coord
        DXw = (DLTproj(pos + dx,coeffs) - DLTproj(pos - dx,coeffs))/(2*norm(dx)); 
        DYw = (DLTproj(pos + dy,coeffs) - DLTproj(pos - dy,coeffs))/(2*norm(dy)); 
        DZw = (DLTproj(pos + dz,coeffs) - DLTproj(pos - dz,coeffs))/(2*norm(dz)); 

        % indices
        x_ind = 2*(ii-1) + 1;
        y_ind = x_ind + 1;

        % assign output
        C(x_ind,1) = DXw(1);
        C(x_ind,2) = DYw(1);
        C(x_ind,3) = DZw(1);
        
        % y coord
        C(y_ind,1) = DXw(2);
        C(y_ind,2) = DYw(2);
        C(y_ind,3) = DZw(2);

    end

end