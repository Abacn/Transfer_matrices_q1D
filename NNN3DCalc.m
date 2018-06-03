% Run Script
% option. 0: calculate xi and rho; 1, only xi; 2, only rho,
% > 4, xi1~option.

function [V, D] = NNN3DCalc(option)

global gcount;
gcount = 0;

if nargin<1
    option = 0;
end

param = load('param3D.dat');
nr = param(1);
nt = param(2);
ns = 100;
dmax = param(3);
betaF = param(4);

if 0==option || 1==option
    fprintf('Calculate eigenvalue\n');
    [V, D] = NNN3DEigs(nr, nt, ns, dmax, betaF);

    dlmwrite('D.dat', D, 'delimiter', '\t', 'Precision', 10);
    % File too large. Nearly 1G for 100x100x100.
    %dlmwrite('V.dat', V, 'delimiter', '\t', 'Precision', 10);

    xi(1) = 1.0 / log(abs(D(1, 1)) / abs(D(2, 2)));
    xi(2) = 1.0 / log(abs(D(1, 1)) / abs(D(3, 3)));
    dlmwrite('xi.dat', xi', 'delimiter', '\t', 'Precision', 10);

end

if any(option==[0,2])  % calculate density
    if 0==option
        rho = NNN3DDensity(nr, nt, ns, dmax, betaF, V, D);
    else
        rho = NNN3DDensity(nr, nt, ns, dmax, betaF);
        V = rho;
    end
    dlmwrite('rho.dat', rho, 'delimiter', '\t', 'Precision', 10);
end

if option >= 4
    vsize = NNN3DMatrix(0, nr, nt, ns, dmax, betaF, 'size');
    D = eigs(@(x) NNN3DMatrix(x, nr, nt, ns, dmax, betaF), vsize, option+1, 'lm');
    [diaD2, ind2] = sort(abs(D),'descend');
    D = diag(diaD2);
    Dsorted = sort(abs(D),'descend');
    for rp=2:option+1
      xi(rp-1) = 1.0 / log(abs(D(1)) / abs(D(rp)));
    end
    dlmwrite('D.dat', diag(D), 'delimiter', '\t', 'Precision', 10);
    dlmwrite('xi.dat', xi', 'delimiter', '\t', 'Precision', 10);
end
if nargout<1
    clear V D;
end
end

% [V,D]=NNN2DCalc('m');