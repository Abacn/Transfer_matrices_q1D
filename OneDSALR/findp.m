% Find pressure, given temperature and density
% Secant method
function p = findp(rho, T, coeffs, tol)
    % start from hard sphere
    if(nargin < 4)
        tol  = 1e-6;
    end
    N = 200;
    sigma = coeffs(1);
    beta = 1/T;
    p1 = rho*T/(1-sigma*rho);
    rho1 = findrho(p1, beta, coeffs, N);
    % First search a border
    if(rho1 < rho)
            p2 = 2*p1;
    elseif(rho1 > rho)
        p2 = p1/2;
    else
        p = p1;
        return;
    end
    rho2 = findrho(p2, beta, coeffs, N);
    subp = p2-p1;
    err = abs(subp/(p2+p1)*2);
    while(err > tol)
        p0 = p1; p1 = p2;
        rho0 = rho1; rho1 = rho2;
        p2 = p1 - (rho1-rho)/(rho1-rho0)*subp;
        if(p2<=0)
            p2 = p1 / 10;
        end
        rho2 = findrho(p2, beta, coeffs, N);
        subp = p2-p1;
        err = abs(subp/(p2+p1)*2);
    end
    p = p2;
end

