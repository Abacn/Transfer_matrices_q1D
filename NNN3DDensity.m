function rho = NNN3DDensity(nr, nt, ns, dmax, betaF, V, D)
    global gcount;
    % flag for eigenvalue and eigenvector input
    egflg = false;
    if nargin==2
      qM = nt(:,1);
      eigM = ns(1,1);
      egflg = true;
    elseif nargin==7
      qM = V(:,1);
      eigM = D(1,1);
      egflg = true;
    end
    if nargin<5 
      param = load('param3D.dat');
        nr = param(1);
        nt = param(2);
        ns = 100;
        dmax = param(3);
        betaF = param(4);
    end
    
    reloadflag =false;
    stofname = 'NNN3DCor_redo';
    
    if exist(strcat(stofname,'.mat'), 'file')
        load(stofname);
        if param_sv(1)==nr && param_sv(2)==nt && param_sv(3)==dmax && param_sv(4)==betaF
            reloadflag = true;
            eigMT = eigM;
        end 
    end
    lenX = NNN3DMatrix(0, nr, nt, ns, dmax, betaF, 'size');
    fprintf('Reduced vector size: %d/%d\n', lenX, nr*nr*nt*ns);
    if (reloadflag || egflg) && lenX ~= size(qM,1)
        warning('Input eigenvector size does not match. Redo calculation.');
        reloadflag = false;
        egflg = false;
    end
    % No reloadable data available, do the eigenvalue calculation.
    if ~reloadflag
      gcount = 0;
      lenX = NNN3DMatrix(0, nr, nt, ns, dmax, betaF, 'size');
      if ~egflg
          fprintf('Calculating q\n');
          [qM, eigM] = eigs(@(x) NNN3DMatrix(x, nr, nt, ns, dmax, betaF), lenX, 1, 'lm');
      end
      gcount = 0;
      fprintf('Calculating q^{-1}\n');
      [qMit, eigMT] = eigs(@(x) NNN3DTMatrix(x, nr, nt, ns, dmax, betaF), lenX, 1, 'lm');

      param_sv = [nr, nt, dmax, betaF];
      save(stofname, 'qM', 'qMit', 'eigM', 'param_sv');
    end
    
    rn = qMit.'*qM;   % reverse normalization
    rho = rn/((qMit .* gentangent(nr, nt, ns, dmax, betaF)).'* qM);
end

% return d K(\Delta x) /d (\beta F), results.*K gives the tangent matrix
function results = gentangent(nr, nt, ns, dmax, betaF)
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
        sigmatmp = 1-rh0-rh1+2*sqrt(rh0*rh1)*ctlist;
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

    results = zeros(totalNS, 1);
    
    jstart = 1;
    for rrh0=1:nr
      r0 = rlist(rrh0);
      for rrh1=1:nr
        for rt0=1:nt
          t0 = tlist(rt0);
          ns1=NStable(rt0, rrh0, rrh1);
          xlist = slist(1:ns1)+sigmatable(rt0, rrh0, rrh1);
          xlist(end) = xlist(end)+1/betaF;
          results(jstart:jstart+ns1-1,:) = xlist.';
          jstart = jstart + ns1;
        end
      end
    end

end