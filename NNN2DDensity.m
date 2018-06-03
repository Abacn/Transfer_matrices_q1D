function rho = NNN2DDensity(nt, ns, ymax, betaF, tanhc, V, D)
    % Called by script, output to the working dir
    % two mode
    % (1) parameter input by argin: NNN2DDensity(nt, ns, ymax, betaF, tanhc, (optional: V, D))
    % (2) parameter input by file param.dat: NNN2DDensity((optional: V, D))
    % if V, D is assigned, used it to get qM and eigM, save half of time.
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
      load param.dat;
      nt = param(1);
      ns = 100;
      ymax = param(2);
      betaF = param(3);
      tanhc = param(4);
    end
    
    reloadflag =false;
    stofname = 'NNN2DCor_redo';
    
    if exist(strcat(stofname,'.mat'), 'file')
        load(stofname);
        if param_sv(1)==nt && param_sv(2)==ymax && param_sv(3)==betaF && param_sv(4)==tanhc
            reloadflag = true;
            eigMT = eigM;
        end 
    end
    
    gcount = 0;
    vsize = NNN2DMatrix(0, nt, ns, ymax, betaF, tanhc, 'size');
    fprintf('Reduced vector size: %d/%d\n', vsize, nt*nt*nt);
    if (reloadflag || egflg) && vsize ~= size(qM,1)
        warning('Input eigenvector size does not match. Redo calculation.');
        egflg = false;
        reloadflag = false;
    end
    if ~(egflg || reloadflag)
      fprintf('Calculating q\n');
      [qM, eigM] = eigs(@(x) NNN2DMatrix(x, nt, ns, ymax, betaF, tanhc), vsize, 1, 'lm');
    end
    gcount = 0;
    if ~reloadflag
      fprintf('Calculating q^{-1}\n');
      [qMit, eigMT] = eigs(@(x) NNN2DTMatrix(x, nt, ns, ymax, betaF, tanhc), vsize, 1, 'lm');
      param_sv = [nt, ymax, betaF, tanhc];
      save(stofname, 'qM', 'qMit', 'eigM', 'param_sv');
    end
    rn = qMit'*qM;   % reverse normalization
    rho = rn/((qMit .* gentangent(nt, ns, ymax, betaF, tanhc)).'* qM);
    
    dlmwrite('rho.dat', rho, 'delimiter', '\t', 'Precision', 10);
    
end
% rho = NNN2DDensity(100, 100, 0.5, 5, 3);


% return d K(\Delta x) /d (\beta F), results.*K gives the tangent matrix
function results = gentangent(nt, ns, ymax, betaF, c)
  deltat = 2.0 / (nt-1);

  if ymax < 1.0
    xmin = sqrt(1.0 - ymax^2);
  else
    xmin = 0.0;
  end
  smaxmax = 1.0 - xmin;
  deltas = smaxmax / (ns-1);
  b = ymax / 2.0;
  a = ymax / 2.0 - b * tanh(c);
  
  % two end all considered. May test if midpoint rule results sth. different.
  tlist = linspace(-1, 1, nt);
  slist = linspace(0, smaxmax, ns); % scaled s. The first element denotes smax, and it is ideal
  ylist = tlist*a + tanh(tlist*c)*b;
  tmp1 = repmat(ylist,nt,1); tmp2 = tmp1-tmp1.'; tmp3 = tmp2.*tmp2; tmp3(tmp3>1.0)=1.0;
  sigmatable = sqrt(1.0-tmp3);
  % the number of s we have
  NStable = ceil((1.0-sigmatable)/deltas)+1;
  NStable(NStable > ns) = ns;
  % total vector length
  totalNS = sum(NStable(:));
  
  results = zeros(totalNS, 1);
  
  istart = 1;
  % three for block stacked
  for rt0=1:nt
      % y0 = ylist(rt0);
      % t0 = tlist(rt0);
      % sechct = sech(c*t0);
      % dydt = a + b*c*sechct*sechct; % use this
      for rt1=1:nt
        % y1= ylist(rt1); % y1 not used. Only for debug usage.
        sigma01 = sigmatable(rt0, rt1); % also multiply this in entries
        ns0 = NStable(rt0, rt1);
        xlist = sigma01+slist(1:ns0);
        xlist(end) = xlist(end)+1/betaF;
        results(istart:istart+ns0-1) = xlist;
        istart = istart + ns0;
      end
  end
end
