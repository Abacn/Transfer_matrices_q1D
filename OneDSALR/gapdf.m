% Calculate the cluster size distribution function
% P - pressure, beta - temperature,
% coeffs = [\sigma, \lambda, \kappa, \epsilon, \xu]
% which defines SALR potential u = \begin{case}
% \infty, (r < \sigma) \\
% -\epsilon, (\sigma <= r < \lambda*\sigma) \\
% \xi\epsilon(\kappa-r/\sigma), (\lambda*\sigma <= r < \kappa*\sigma) \\
% 0, (r > \kappa*\sigma)
% set 2*pi*m/h^2 = 1
% consider third nearest neighbor
% N: up to cluster of size N
% Outputs:
% rlist: a list of particle distances s; pros: probability density of the neighbor at distance s


function [rlist, pros] = gapdf(P, beta, coeffs, divide)
    if( nargin<4)
        divide = 200;
    end
    divide = divide + 1; % tail
    sigma = coeffs(1);
    kappasigma = coeffs(3)*sigma;
    lambdasigma = coeffs(2)*sigma;
    ds = (kappasigma-sigma)/(divide-1);
    % divide2 = floor((kappasigma-sigma)/ds);
    rlist = linspace(sigma, kappasigma, divide).'+ds/2; % Magic. Increased the accurancy pretty much
    pot = potential(rlist, coeffs);
    pot = pot + P*(rlist-sigma);
    tmat = pot*ones(1,divide);
    % tmat: s(i) by s(i+1)
    for rp=1:divide
        % plus next nearest neighbor
        tmat(:,rp) = exp(-beta*(tmat(:,rp) + potential(rlist+rlist(rp), ...
            coeffs)));
    end
    bpks = beta*P*(coeffs(3)-1)*sigma;  % contribution of tail
    tail = exp(-bpks)/(beta*P);
    tmat(divide, 1:divide) = tail/ds; % add tail to the last row
    mlen = divide*divide;
    % Generating 3NN transfer matrix t2mat (sparse) according to
    % NNN transfer matrix M
    vlen = mlen*divide;  % vector length
    xlist = zeros(vlen,1);ylist = zeros(vlen,1);elist = zeros(vlen,1);
    M2 = zeros(divide, mlen);
    bind = 1:divide;
    for rp=bind
        for rq=bind
            xindex_head = (rq-1)*divide;  % first x index
            yindex = (rp-1)*divide+rq;    % y indeces
            lindex = (yindex-1)*divide; % first store index x/y/e in array
            xlist(bind+lindex) = bind+xindex_head;
            ylist(bind+lindex) = yindex;
            values = tmat(:,rq).*exp(-beta*potential(rlist+rlist(rq)+rlist(rp), coeffs));
            elist(bind+lindex) = values;
            M2(:,yindex) = values;
        end
    end
    M = sparse(xlist, ylist, elist, mlen, mlen);
    [qM, eigM] = eigs(M,1);
    [qMit, ~] = eigs(M.',1);
    rn = qMit'*qM*eigM;   % reverse normalization, denominator
    qMitmat = reshape(qMit, [divide, divide]);
    for rp=bind
        xindex_head = (rp-1)*divide;
        M2(:, xindex_head+bind) = qMitmat .* M2(:, xindex_head+bind);
    end
    pros = M2 * qM / (rn*ds);
    if(nargout<1)
        plot(rlist(1:end-1), pros(1:end-1));
        rho = findrho(P, beta, coeffs);
        rlist = rho;
        pros = [];
    end
end
