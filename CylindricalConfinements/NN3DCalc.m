function [xi, rho, V, D] = NN3DCalc(param)
    if nargin<1
        % test
        param = load('param3D.dat');
    end
    nr = param(1);
    nt = param(2);
    dmax = param(3);
    betaF = param(4);

    [V, D, tmat] = NN3DEigs(nr, nt, dmax, betaF);
    % calculate xi
    xi(1) = 1.0 / log(abs(D(1, 1)) / abs(D(2, 2)));
    xi(2) = 1.0 / log(abs(D(1, 1)) / abs(D(3, 3)));
    % calculate rho
    eigM = D(1,1);
    qM = V(:,1);
    [qMit, eigMit] = eigs(tmat.', 1);
    
    rho = qMit.'*qM*eigM / (qMit.'*NNN2DDifMat(nr, nt, dmax, betaF)*qM);
    
    dlmwrite('xi.dat', xi', 'delimiter', '\t', 'Precision', 10);
    dlmwrite('rho.dat', rho, 'delimiter', '\t', 'Precision', 10);
    
    dlmwrite('D.dat', D, 'delimiter', '\t', 'Precision', 10);
    dlmwrite('V.dat', V, 'delimiter', '\t', 'Precision', 10);
    if nargout<1
        clear('xi');
    end
end

%dt=0.5;p=4;[xi, rho] = NN3DCalc([100,100,sqrt(dt),p]);

function DifMat = NNN2DDifMat(nr, nt, dmax, betaF)

% The particle diameter is set to be 1.0.
% Nearest neighbor transfer matrix eigenvalues

% ***  2*drh*dt/2/betaF is not included in the matrix ***

if dmax>sqrt(3)/2
    error('NNN interaction possible.');
end

rhmax = (dmax/2)^2;
deltarh = rhmax/nr;
rhlist = linspace(rhmax-deltarh, 0, nr) + deltarh/2; % midpoint rule
deltat = pi/nt;
tlist = linspace(0, pi-deltat, nt);

dttable = repmat(linspace(1,nt,nt).', 1, nt);
dtexpanded = dttable(:);

DifMat = zeros(nr*nt, nr*nt);

for rp=1:nr
  xstart = nt*(rp-1)+1;
  for rq=rp:nr
    ystart = nt*(rq-1)+1;
    rh0 = rhlist(rp);
    rh1 = rhlist(rq);
    sigmalist = sqrt(1-rh0-rh1+2*sqrt(rh0*rh1)*cos(tlist));
    blist = exp(-betaF*sigmalist) .* sigmalist .* (1+1./(betaF*sigmalist));
    bmat = reshape(blist(dtexpanded), nt, nt);
    DifMat(xstart:xstart+nt-1, ystart:ystart+nt-1) = bmat;
    if rp ~= rq
      DifMat(ystart:ystart+nt-1, xstart:xstart+nt-1) = bmat;
    end
  end
end

end
