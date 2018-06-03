function gs = NN3DCor(uplimit, param)
    global gcount;
    option = 'xy';
    %if nargin<2, option = 'xy';end
    if nargin<1 || uplimit<6, uplimit = 20;end % uplimit at least be 5
    nxis = 2;  % number of correlations investigated
    gcount = 0;

    if nargin<2
        % test
        param = load('param3D.dat');
    end
    nr = param(1);
    nt = param(2);
    dmax = param(3);
    betaF = param(4);
    
    ysxs = genvector(nr, nt, dmax);
    lenX = size(ysxs, 1);
    
    % try to reload results, if any
    reloadflag = false;
    stofname = 'NN3DCor_redo';
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
        fprintf('Calculating eigenvectors\n');
        [V, D, tmat] = NN3DEigs(nr, nt, dmax, betaF);
        qM = V(:,1);
        eigM = D(1,1);
        [qMit, eigMT]=eigs(tmat.', 1);
        param_sv = param;
        save(stofname, 'qM', 'qMit', 'eigM', 'tmat', 'param_sv');
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
            rightv(:, rq) = tmat * rightv(:, rq) /eigM;
            gs(rp, rq) = leftv(rq, :) * rightv(:, rq) - gs0(rq);
        end
    end
    
    gs = [gs0; gs];
    % write files
    % gs.dat outputs <x_i>^2 and <x_i x_(i+m)> - <x_i>^2 list, m from 0 to uplimit
    dlmwrite('gs.dat', gs, 'delimiter', '\t', 'Precision', 10); 
end

function results = genvector(nr, nt, dmax)
    rhmax = (dmax/2)^2;
    deltarh = rhmax/nr;
    rhlist = linspace(rhmax-deltarh, 0, nr) + deltarh/2; % midpoint rule
    rlist = sqrt(rhlist);
    deltat = pi/nt;
    tlist = linspace(0, pi-deltat, nt);

    results = zeros(nt*nr, 2);
    
    % (r, theta)
    xstart = 1;
    for rp=1:nr
      for rq=rp:nr
        results(xstart, :) = [rlist(rp), tlist(rq)];
        xstart = xstart + 1;
      end
    end

end