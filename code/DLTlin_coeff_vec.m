function [vec_x,vec_y] = DLTlin_coeff_vec(pos,coeffs)
% compute vector of mapping function gradients wrt to coefficients
% returns 2 column vectors of size ncoeff x 1

% set up problem
ncoeff = size(coeffs,1);
vec_x = NaN(ncoeff,1);
vec_y = NaN(ncoeff,1);
pos_mod = [pos 1];

% pre-compute denominator numerator
denom = coeffs(9)*pos(1) + coeffs(10)*pos(2) + coeffs(11)*pos(3) + 1;
num_x = coeffs(1)*pos(1) + coeffs(2)*pos(2) + coeffs(3)*pos(3) + coeffs(4);
num_y = coeffs(5)*pos(1) + coeffs(6)*pos(2) + coeffs(7)*pos(3) + coeffs(8);
denom_sq = denom^2;

% loop over coefficients
    for ii = 1:ncoeff
        
        if 1 <= ii && ii <= 4
            vec_y(ii) = 0;
            vec_x(ii) = pos_mod(ii)/denom;
        elseif 5 <= ii && ii <= 8
            vec_x(ii) = 0;
            vec_y(ii) = (denom*pos_mod(ii-4) - num_y*0)/denom_sq;
        elseif 9 <= ii && ii <= 11
            vec_x(ii) = (denom*0 - num_x*pos(ii-8))/denom_sq;
            vec_y(ii) = (denom*0 - num_y*pos(ii-8))/denom_sq;
        end

    end

end
