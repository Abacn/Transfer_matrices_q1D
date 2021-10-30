% Calculate the partition function (\Zeta/N)
% P - pressure, beta - temperature,
% coeffs = [\sigma, \lambda, \kappa, \epsilon, \xi]
% which defines SALR potential u = \begin{case}
% \infty, (r < \sigma) \\
% -\epsilon, (\sigma <= r < \lambda*\sigma) \\
% \xi\epsilon(\kappa-r/\sigma), (\lambda*\sigma <= r < \kappa*\sigma) \\
% 0, (r > \kappa*\sigma)
% set 2*pi*m/h^2 = 1
% consider third next nearest neighbor
% Use two part discretization: half of rlist is 1~??1 , another half > ??1
function zeta = Pfunc_isobaric_3NN_dv(P, beta, coeffs, divide)
% set 2*pi*m/h^2 = 1
if ( nargin<4)
    divide = 300;
end
divide = divide + 1; % add tail
sigma = coeffs(1);
lambdasigma = coeffs(2)*sigma;
kappasigma = coeffs(3)*sigma;

optdflag = false;        % optimizing discretization flag
dPoint = lambdasigma-sigma;
% The margin (+0) can change. Ex: +1/8, result in different oscillation
% period. If it is 0, period->infinity and the result nearly monotonically
% reach to the true value
if dPoint > sigma
    divideA = floor(divide/2);
    divideB = divide - divideA;
    divide = divideA + divideB;
    % let error from the discontinuity always cancelled
    dsA = (dPoint-sigma)/(divideA+0.5);
    RealdPoint = dPoint - dsA/2;
    dsB = (kappasigma-RealdPoint)/(divideB-1);
    dsBA = dsB/dsA;
    dsBAlist = [ones([divideA,1]); ones([divideB-1, 1]).*dsBA; 1];
    rlistA = linspace(sigma, RealdPoint-dsA, divideA).'+dsA/2;
    rlistB = linspace(RealdPoint, kappasigma, divideB).' + dsB/2;
    rlist = [rlistA; rlistB];
    optdflag = true;
    ds = dsA;
else
    ds = (kappasigma-sigma)/(divide-1);
    rlist = linspace(sigma, kappasigma, divide).'+ds*1/2; % +ds/2 % Magic. Increased the accurancy pretty much
end
pot = potential(rlist, coeffs);
pot = pot + P*(rlist-sigma);
tmat = pot*ones(1,divide);
% tmat: s(i) by s(i+1)
for rp=1:divide
% plus next nearest neighbor
tmat(:,rp) = exp(-beta*(tmat(:,rp) + potential(rlist+rlist(rp), ...
                                               coeffs)));
end
if optdflag
    tmat = tmat .* (dsBAlist * ones(1, length(dsBAlist)));
end
bpks = beta*P*(coeffs(3)-1)*sigma;  % contribution of tail
tail = exp(-bpks)/(beta*P);
tmat(divide, 1:divide) = tail/ds; % add tail to the last row
mlen = divide*divide;

% Generating 3NN transfer matrix t2mat (sparse) according to
% NNN transfer matrix tmat
vlen = mlen*divide;  % vector length
xlist = zeros(vlen,1);
ylist = zeros(vlen,1);
elist = zeros(vlen,1);
bind = 1:divide;
for rp=bind
for rq=bind
xindex_head = (rq-1)*divide;  % first x index
yindex = (rp-1)*divide+rq;    % y indeces
lindex = (yindex-1)*divide; % first store index x/y/e in array
xlist(bind+lindex) = bind+xindex_head;
ylist(bind+lindex) = yindex;
elist(bind+lindex) = tmat(:,rq) ...
.*exp(-beta*potential(rlist+rlist(rq)+rlist(rp), coeffs));
end
end
t2mat = sparse(xlist, ylist, elist, mlen, mlen);
lmax = eigs(t2mat,1);
zeta = sqrt(1/beta)*lmax*ds*exp(-beta*P*sigma);
end

