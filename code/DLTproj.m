function Xc = DLTproj(pos,coeffs)
% project position onto camera view

    Xc = NaN(1,2);
    denom = coeffs(9)*pos(1) + coeffs(10)*pos(2) + coeffs(11)*pos(3) + 1;
    Xc(1) = (coeffs(1)*pos(1) + coeffs(2)*pos(2) + coeffs(3)*pos(3) + coeffs(4))/denom;
    Xc(2) = (coeffs(5)*pos(1) + coeffs(6)*pos(2) + coeffs(7)*pos(3) + coeffs(8))/denom;

end

