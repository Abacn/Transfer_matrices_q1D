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
% Use simpson sum to construct transfer matrix
function zeta = Pfunc_isobaric_3NN_sp(P, beta, coeffs, divide)
% Simpson rule for numerical integration
% set 2*pi*m/h^2 = 1
if ( nargin<4)
    divide = 300;
end
divideB = 2*divide + 1; % add tail
sigma = coeffs(1);
lambdasigma = coeffs(2)*sigma;
kappasigma = coeffs(3)*sigma;

ds = (kappasigma-sigma)/divide;
dsB = (kappasigma-sigma)/(divideB-1);

% First construct a matrix sizes (3*m+1)*(3*m+1)
rlistB = linspace(sigma, kappasigma, divideB).';
rlist = rlistB(2:2:end);
rlist(end+1) = rlistB(end) + ds;
pot = potential(rlistB, coeffs);
pot = pot + P*(rlistB-sigma);
tmatB = pot*ones(1,divideB);
% tmat: s(i) by s(i+1)
for rp=1:divideB
% plus next nearest neighbor
    tmatB(:,rp) = exp(-beta*(tmatB(:,rp) + potential(rlistB+rlistB(rp), ...
                                            coeffs)));
end
% Now recover the matrix sizes
divide = divide + 1;  % add tail
tmat = zeros(divide, divide);
tmpind = 1:divide-1;
% 2D simpson rule for entries
for rp = tmpind
    tmat(tmpind, rp) = ( ...
        tmatB(tmpind*2-1, rp*2-1) + 4*tmatB(tmpind*2-1, rp*2) + tmatB(tmpind*2-1, rp*2+1) ...
      + 4*tmatB(tmpind*2, rp*2-1) + 16*tmatB(tmpind*2, rp*2) + 4*tmatB(tmpind*2, rp*2+1) ...
      + tmatB(tmpind*2+1, rp*2-1) + 4*tmatB(tmpind*2+1, rp*2) + tmatB(tmpind*2+1, rp*2+1) ...
    )/36;
end
% 1D simpson rule for NNN tail
tmat(tmpind, divide) = (tmatB(tmpind*2-1, divideB) + 4*tmatB(tmpind*2, divideB) + tmatB(tmpind*2+1, divideB))/6;

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

