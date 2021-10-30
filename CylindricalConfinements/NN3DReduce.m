function [xi, V, D, tmat] = NN3DReduce(param)
    if nargin<1
        % test
        param = load('param3D.dat');
    end
    nr = param(1);
    nt = param(2);
    dmax = param(3);
    betaF = param(4);

    [V, D, tmat] = NN3DReduceMatrix(nr, nt, dmax, betaF);
    % calculate xi
    xi(1) = 1.0 / log(abs(D(1, 1)) / abs(D(2, 2)));
    xi(2) = 1.0 / log(abs(D(1, 1)) / abs(D(3, 3)));
    % calculate rho
    eigM = D(1,1);
    qM = V(:,1);
end

function [V, D, tmat] = NN3DReduceMatrix(nr, nt, dmax, betaF)

% Angular part is reduced

% ***  2*drh*dt/2/betaF not included in the matrix ***

if dmax>sqrt(3)/2
    error('NNN interaction possible.');
end

rhmax = (dmax/2)^2;
deltarh = rhmax/nr;
rhlist = linspace(rhmax-deltarh, 0, nr) + deltarh/2; % midpoint rule
deltat = pi/nt;
tlist = linspace(0, pi-deltat, nt);

tmat = zeros(nr, nr);

for rp=1:nr
  xstart = nt*(rp-1)+1;
  for rq=rp:nr
    ystart = nt*(rq-1)+1;
    rh0 = rhlist(rp);
    rh1 = rhlist(rq);
    
    sigmalist = sqrt(1-rh0-rh1+2*sqrt(rh0*rh1)*cos(tlist));
    blist = exp(-betaF*sigmalist);
    
    bsum=sum(blist);
    tmat(rp, rq) = bsum;
    if rp ~= rq
      tmat(rq, rp) = bsum;
    end
  end
end


[V, D]=eigs(tmat, 3);
diaD = diag(D);
[diaD2, ind2] = sort(abs(diaD),'descend');
D = diag(diaD2);
V = V(:,ind2);
%prefx = 0.5*deltarh*deltat/betaF;
%D = D*prefix;
end


% [V,D]=NN3DEigs(50,100,0.5,5);