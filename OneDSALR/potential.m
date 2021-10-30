% SALR potential
function pot = potential(rs, coeffs)
    % coeffs = [\sigma, \lambda, \kappa, \epsilon, \zeta]
    sigma = coeffs(1);
    pot = zeros(size(rs));
    index = 1;
    for r=abs(rs(:)).'
        if(r<sigma) 
            u = inf;
        elseif (r<coeffs(2)*sigma)
            u = -coeffs(4);
        elseif (r<coeffs(3)*sigma)
            u = coeffs(5)*coeffs(4)*(coeffs(3)-r/sigma);
        else
            u = 0;
        end
        pot(index) = u;
        index = index + 1;
    end
end
