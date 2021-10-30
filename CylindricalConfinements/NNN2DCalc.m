% option: 0-calculate first two xi.
% >=4, the number of xi that will be calculated

function [V, D] = NNN2DCalc(option)

if nargin < 1
    option = 0;
end
global gcount;
gcount = 0;

load param.dat;
nt = param(1);
ns = 100;
ymax = param(2);
betaF = param(3);
tanhc = param(4);

if 0==option
	[V, D] = NNN2DEigs(nt, ns, ymax, betaF, tanhc);

	dlmwrite('D.dat', D, 'delimiter', '\t', 'Precision', 10);
	dlmwrite('V.dat', V, 'delimiter', '\t', 'Precision', 10);

	xi(1) = 1.0 / log(abs(D(1, 1)) / abs(D(2, 2)));
	xi(2) = 1.0 / log(abs(D(1, 1)) / abs(D(3, 3)));
	dlmwrite('xi.dat', xi', 'delimiter', '\t', 'Precision', 10);
elseif option >= 4
	vsize = NNN2DMatrix(0, nt, ns, ymax, betaF, tanhc, 'size');
    D = eigs(@(x) NNN2DMatrix(x, nt, ns, ymax, betaF, tanhc), vsize, option+1, 'lm');
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