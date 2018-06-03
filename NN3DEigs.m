function [V, D, tmat] = NN3DEigs(nr, nt, dmax, betaF)

% The particle diameter is set to be 1.0.
% Nearest neighbor transfer matrix eigenvalues

% ***  2*drh*dt/2/betaF not included in the matrix ***

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

tmat = zeros(nr*nt, nr*nt);

for rp=1:nr
  xstart = nt*(rp-1)+1;
  for rq=rp:nr
    ystart = nt*(rq-1)+1;
    rh0 = rhlist(rp);
    rh1 = rhlist(rq);
    sigmalist = sqrt(1-rh0-rh1+2*sqrt(rh0*rh1)*cos(tlist));
    blist = exp(-betaF*sigmalist);
    bmat = reshape(blist(dtexpanded), nt, nt);
    tmat(xstart:xstart+nt-1, ystart:ystart+nt-1) = bmat;
    if rp ~= rq
      tmat(ystart:ystart+nt-1, xstart:xstart+nt-1) = bmat;
    end
  end
end

[V, D]=eigs(tmat, 5);
diaD = diag(D);
[diaD2, ind2] = sort(abs(diaD),'descend');
D = diag(diaD(ind2));
V = V(:,ind2);
%prefx = 0.5*deltarh*deltat/betaF;
%D = D*prefix;
end


% [V,D]=NN3DEigs(50,100,0.5,5);