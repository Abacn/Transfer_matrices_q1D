function [V, D] = NNN2DEigs(nt, ns, ymax, betaF, tanhc, option)

global gcount;
% The particle diameter is set to be 1.0.
% Nearest neighbor transfer matrix eigenvalues

% *** dy = ymax/(ny-1); ***
if nargin<6
  option = 'm';
end
if ymax>1.5
    error('3NN interaction possible.');
end
deltat = 2.0 / (nt-1);
if ymax < 1.0
    xmin = sqrt(1.0 - ymax^2);
else
    xmin = 0.0;
end
smaxmax = 1.0 - xmin;
if strcmpi(option, 'm')
    gcount = 0;
    vsize = NNN2DMatrix(0, nt, ns, ymax, betaF, tanhc, 'size');
    fprintf('Reduced vector size: %d/%d\n', vsize, nt*nt*nt);
    [V, D] = eigs(@(x) NNN2DMatrix(x, nt, ns, ymax, betaF, tanhc), vsize, 3, 'lm');
    % D = D*exp(-betaF*xmin);
else
    error('Unknown Option "%s"', option);
end
% prefx = deltat*deltas;
% D = D*prefx;
% adjust eigenvalue by magnitude
diaD = diag(D);
[diaD2, ind2] = sort(abs(diaD),'descend');
D = diag(diaD(ind2));
V = V(:,ind2);
end

% Nx =  NNN2DMatrix(0, 100, 100, 0.7, 7.0, 3, 'size'); X = ones(Nx,1);Y = NNN2DMatrix(X, 100, 100, 0.7, 7.0, 3);
% lambda1 = eigs(@(X) NNN2DMatrix(X, 100, 100, 0.7, 7.0, 3), Nx, 3, 'LM')