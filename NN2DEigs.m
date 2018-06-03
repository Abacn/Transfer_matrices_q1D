function [V, D] = NN2DEigs(nt, ymax, betaF, c)

% The particle diameter is set to be 1.0.
% Nearest neighbor transfer matrix eigenvalues

% *** dy = ymax/(ny-1); ***

if ymax>sqrt(3)/2
    error('NNN interaction possible.');
end
b = ymax / 2.0;
a = ymax / 2.0 - b * tanh(c);

deltat = 2.0 / (nt-1);
tlist = linspace(-1, 1, nt);
ylist = tlist*a + tanh(tlist*c)*b;
tmp1 = repmat(ylist,nt,1); tmp2 = tmp1-tmp1.'; tmp3 = tmp2.*tmp2; tmp3(tmp3>1.0)=1.0;
sigmatable = sqrt(1.0-tmp3);
xmin = sqrt(1.0 - ymax);
sechct = sech(c*tlist);
dydt = a + b*c*sechct.^2; % dy/dt
tmat = exp(-betaF*(sigmatable-xmin)) .* repmat(dydt, size(sigmatable, 1), 1);
[V, D]=eigs(tmat, 5);
prefx = exp(-betaF*xmin)*deltat;
D = D*prefx;
end

% Nx =  NNN2DMatrix(0, 100, 100, 0.7, 7.0, 3, 'size'); X = ones(Nx,1);Y = NNN2DMatrix(X, 100, 100, 0.7, 7.0, 3);
% lambda1 = eigs(@(X) NNN2DMatrix(X, 100, 100, 0.7, 7.0, 3), Nx, 3, 'LM')
