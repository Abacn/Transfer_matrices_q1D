function [V, D] = NNN3DEigs(nr, nt, ns, dmax, betaF)
% The particle diameter is set to be 1.0.
% Next Nearest neighbor transfer matrix eigenvalues

% ***  drh*dt/2/betaF not included in the matrix ***
global gcount;
if dmax>1
    error('3NN interaction possible.');
end

% first estimate the vector size
vsize = NNN3DMatrix(0, nr, nt, ns, dmax, betaF, 'size');
fprintf('Reduced vector size: %d/%d\n', vsize, nr*nr*nt*ns);
gcount = 0;

[V, D] = eigs(@(x) NNN3DMatrix(x, nr, nt, ns, dmax, betaF), vsize, 3, 'lm');
diaD = diag(D);
[diaD2, ind2] = sort(abs(diaD),'descend');
D = diag(diaD(ind2));
V = V(:,ind2);

end