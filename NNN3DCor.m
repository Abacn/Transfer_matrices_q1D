function gs = NNN3DCor(uplimit, param)
    global gcount;
    option = 'xy';
    %if nargin<2, option = 'xy';end
    if nargin<1 || uplimit<6, uplimit = 20;end % uplimit at least be 5
    nxis = 3;  % number of correlations investigated
    gcount = 0;

    if nargin<2
        % test
        param = load('param3D.dat');
    end
    nr = param(1);
    nt = param(2);
    dmax = param(3);
    betaF = param(4);
    ns = 100;
    ysxs = genvector(nr, nt, ns, dmax, betaF);
    lenX = size(ysxs, 1);
    
    % try to reload results, if any
    reloadflag = false;
    stofname = 'NNN3DCor_redo';
    if exist(strcat(stofname,'.mat'), 'file')
        load(stofname);
        if param_sv(1)==nr && param_sv(2)==nt && param_sv(3)==dmax && param_sv(4)==betaF
            reloadflag = true;
            eigMT = eigM;
        end 
    end
    % No reloadable data available, do the eigenvalue calculation.
    if ~reloadflag
        gcount = 0;
        fprintf('Calculating q\n');
        [qM, eigM] = eigs(@(x) NNN3DMatrix(x, nr, nt, ns, dmax, betaF), lenX, 1, 'lm');
        gcount = 0;
        fprintf('Calculating q^{-1}\n');
        [qMit, eigMT] = eigs(@(x) NNN3DTMatrix(x, nr, nt, ns, dmax, betaF), lenX, 1, 'lm');
        param_sv = param;
        save(stofname, 'qM', 'qMit', 'eigM', 'param_sv');
    end
    gcount = 0;
    fprintf('Calculating correlation\n');
    qMit = qMit.';
    ysxsit = ysxs.';
    denom = qMit*qM;  % denominator, used to normalize qM
    qM = qM/denom;
    leftv = ysxsit .* repmat(qMit, nxis, 1);
    rightv = ysxs .* repmat(qM, 1, nxis);
    gs = zeros(uplimit, nxis);
    gs0 = zeros(1, nxis);
    
    % first <x^2> and <x>^2
    for rq=1:nxis
        gs0(rq) = (qMit * rightv(:, rq))^2;
        gs(1, rq) = leftv(rq, :) * rightv(:, rq) - gs0(rq);
    end
    % then calculate <x_i x_j> - <x^2>
    for rp=2:uplimit+1
        for rq=1:nxis
            rightv(:, rq) = NNN3DMatrix(rightv(:, rq), nr, nt, ns, dmax, betaF) /eigM;
            gs(rp, rq) = leftv(rq, :) * rightv(:, rq) - gs0(rq);
        end
    end
    
    gs = [gs0; gs];
    % write files
    % gs.dat outputs <x_i>^2 and <x_i x_(i+m)> - <x_i>^2 list, m from 0 to uplimit
    dlmwrite('gs.dat', gs, 'delimiter', '\t', 'Precision', 10); 
end

% return (r, theta, s list)
function results = genvector(nr, nt, ns, dmax, betaF)
    rhmax = (dmax/2)^2;
    deltarh = rhmax/nr;
    rhlist = linspace(rhmax-deltarh, 0, nr) + deltarh/2; % midpoint rule
    rlist = sqrt(rhlist);
    deltat = 2*pi/nt;
    tlist = linspace(0, 2*pi-deltat, nt);
    tlist(tlist >= pi) = tlist(tlist >= pi) - 2*pi;      % (-pi, pi)
    ctlist = cos(tlist);
    if dmax < 1.0
        xmin = sqrt(1.0 - dmax^2);
    else
        xmin = 0.0;
    end
    smaxmax = 1.0 - xmin;
    deltas = smaxmax / (ns-1);
    slist = linspace(0, smaxmax, ns);

    btable = zeros(nt, nr, nr);      % table of Boltzmann weights
    sigmatable = zeros(nt, nr, nr);  % table of sigma

    for rp=1:nr
      for rq=rp:nr
        rh0 = rhlist(rp);
        rh1 = rhlist(rq);
        sigmatmp = 1-rh0-rh1+2*sqrt(rh0*rh1)*cos(tlist);
        sigmatmp(sigmatmp < 0) = 0;
        sigmalist = sqrt(sigmatmp);
        blist = exp(-betaF*(sigmalist-xmin));
        btable(:, rp, rq) = blist;
        sigmatable(:, rp, rq) = sigmalist;
        if rp ~= rq
          btable(:, rq, rp) = blist;
          sigmatable(:, rq, rp) = sigmalist;
        end
      end
    end
    NStable = ceil((1.0-sigmatable)/deltas)+1;
    NStable(NStable > ns) = ns;
    totalNS = sum(NStable(:));

    results = zeros(totalNS, 3);
    
    jstart = 1;
    for rrh0=1:nr
      r0 = rlist(rrh0);
      for rrh1=1:nr
        for rt0=1:nt
          t0 = tlist(rt0);
          ns1=NStable(rt0, rrh0, rrh1);
          xlist = slist(1:ns1)+sigmatable(rt0, rrh0, rrh1);
          xlist(end) = (1+xlist(end)*betaF)/betaF;
          results(jstart:jstart+ns1-1,:) = [repmat([r0, t0], ns1, 1), xlist.'];
          jstart = jstart + ns1;
        end
      end
    end

end