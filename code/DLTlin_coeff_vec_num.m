function [vec_x,vec_y] = DLTlin_coeff_vec_num(pos,coeffs)
% compute vector of mapping function gradients wrt to coefficients
% returns 2 column vectors of size ncoeff x 1

% set up problem
ncoeff = size(coeffs,1);
vec_x = NaN(ncoeff,1);
vec_y = NaN(ncoeff,1);

% coefficient perturbation
dCoeff = abs(diag(coeffs*0.0001));

% loop over coefficients
    for ii = 1:ncoeff
        
        % perturb
        dC = (DLTproj(pos,coeffs+dCoeff(:,ii)) - DLTproj(pos,coeffs-dCoeff(:,ii)))/(2*norm(dCoeff(:,ii))); 
        
        % assign output
        vec_x(ii) = dC(1);
        vec_y(ii) = dC(2);
        
    end

end
