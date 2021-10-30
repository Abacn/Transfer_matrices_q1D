% Calculate the cluster size distribution function
% P - pressure, beta - temperature,
% coeffs = [\sigma, \lambda, \kappa, \epsilon, \xi]
% which defines SALR potential u = \begin{case}
% \infty, (r < \sigma) \\
% -\epsilon, (\sigma <= r < \lambda*\sigma) \\
% \xi\epsilon(\kappa-r/\sigma), (\lambda*\sigma <= r < \kappa*\sigma) \\
% 0, (r > \kappa*\sigma)
% set 2*pi*m/h^2 = 1
% consider third next nearest neighbor
% N: up to cluster of size N
% Outputs:
% ks: cdf(n); pros: marginal probabilities of finding n neighbored particle
% that belongs to a cluster
% M_lambda

function [ks, pros] = cdf(P, beta, coeffs, divide, N, eigtol)
    if(nargin<6)
        eigtol = 1e-6;
    end
    if(nargin<5)
        N = 10;
        if ( nargin<4)
            divide = 300;
        end
    end
    divide = divide + 1; % tail
    pros = zeros(N+2,1);
    ks = zeros(N,1);
    [M, M2] = genmat(P, beta, coeffs, divide);
    [qM, eigM] = eigs(M,1);
    [qMit, ~] = eigs(M.',1);
    rn = qMit'*qM;   % reverse normalization
    % M = M/eigM;
    M2 = M2/eigM;
    MP = qM;
    eigM2 = eigs(M2,1);
    eigflag = 0;
    pros(1) = 1;
    for rp=2:N+2
        if(eigflag)
            pros(rp) = pros(rp-1)*eigM2;
        else
            MP = M2*MP;
            pros(rp) = qMit'*MP/rn;
            err = pros(rp)/(pros(rp-1)*eigM2)-1;
            if(abs(err)<eigtol)
                eigflag = 1;
            end
        end
    end
    for rp=1:N
        ks(rp) = rp*(pros(rp)+pros(rp+2)-2*pros(rp+1));
    end
end

function [M, M2] = genmat(P, beta, coeffs, divide)
    if (nargin<4)
        divide = 100;
    end
    
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
        divide2 = divideA + round((lambdasigma-dPoint)/dsB);
        optdflag = true;
        ds = dsA;
    else
        ds = (kappasigma-sigma)/(divide-1);
        rlist = linspace(sigma, kappasigma, divide).'+ds*1/2; % +ds/2 % Magic. Increased the accurancy pretty much
        divide2 = round((lambdasigma-sigma)/ds);
    end

%{
    % Old code
    ds = (kappasigma-sigma)/(divide-1);
    divide2 = floor((lambdasigma-sigma)/ds);
    rlist = linspace(sigma, kappasigma, divide).'+ds/2; % Magic. Increased the accurancy pretty much
%}
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
    % NNN transfer matrix M
    vlen = mlen*divide;  % vector length
    vlen2 = mlen*divide2;
    xlist = zeros(vlen,1);ylist = zeros(vlen,1);elist = zeros(vlen,1);
    xlist2 = zeros(vlen2,1);ylist2 = zeros(vlen2,1);elist2 = zeros(vlen2,1);
    bind = 1:divide;
    bind2 = 1:divide2;
    for rp=bind
        for rq=bind
            xindex_head = (rq-1)*divide;  % first x index
            yindex = (rp-1)*divide+rq;    % y indeces
            lindex = (yindex-1)*divide; % first store index x/y/e in array
            xlist(bind+lindex) = bind+xindex_head;
            ylist(bind+lindex) = yindex;
            elist(bind+lindex) = tmat(:,rq) ...
              .*exp(-beta*potential(rlist+rlist(rq)+rlist(rp), coeffs));
            % Cut the element of tmat where s(i)>lambda*sigma
            lindex2 = (yindex-1)*divide2;
            xlist2(bind2+lindex2) = bind2+xindex_head;
            ylist2(bind2+lindex2) = yindex;
            elist2(bind2+lindex2) = elist(lindex+1:lindex+divide2);
        end
    end
    M = sparse(xlist, ylist, elist, mlen, mlen);
    M2 = sparse(xlist2, ylist2, elist2, mlen, mlen);
end
