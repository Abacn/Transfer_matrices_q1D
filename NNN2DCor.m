% 2D correlation function
% input:
% uplimit: compute up to <x_i x_(i+uplimit)>.

function gs = NNN2DCor(uplimit, param)
    global gcount;
    option = 'xy';
    %if nargin<2, option = 'xy';end
    if nargin<1 || uplimit<6, uplimit = 20;end % uplimit at least be 5
    nxis = 2;  % number of correlations investigated
    gcount = 0;
    if nargin<2
      load param.dat;
    end
    nt = param(1);
    ns = 100;
    ymax = param(2);
    betaF = param(3);
    tanhc = param(4);
    % vectors are [ys, xs]
    ysxs = genvector(nt, ns, ymax, betaF, tanhc);
    lenX = size(ysxs, 1);
    
    % try to reload results, if any
    reloadflag = false;
    stofname = 'NNN2DCor_redo';
    if exist(strcat(stofname,'.mat'), 'file')
        load(stofname);
        if param_sv(1)==nt && param_sv(2)==ymax && param_sv(3)==betaF && param_sv(4)==tanhc
            reloadflag = true;
            eigMT = eigM;
        end 
    end
    % No reloadable data available, do the eigenvalue calculation.
    if ~reloadflag
        gcount = 0;
        [qM, eigM] = eigs(@(x) NNN2DMatrix(x, nt, ns, ymax, betaF, tanhc), lenX, 1, 'lm');
        gcount = 0;
        fprintf('Calculating q^{-1}\n');
        [qMit, eigMT] = eigs(@(x) NNN2DTMatrix(x, nt, ns, ymax, betaF, tanhc), lenX, 1, 'lm');
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
            rightv(:, rq) = NNN2DMatrix(rightv(:, rq), nt, ns, ymax, betaF, tanhc)/eigM;
            gs(rp, rq) = leftv(rq, :) * rightv(:, rq) - gs0(rq);
        end
    end
    
    gs = [gs0; gs];
    % write files
    % gs.dat outputs <x_i>^2 and <x_i x_(i+m)> - <x_i>^2 list, m from 0 to uplimit
    dlmwrite('gs.dat', gs, 'delimiter', '\t', 'Precision', 10); 
end

% return [y(y0, y1, s0), \Delta x(y0, y1, s0)] 
function results = genvector(nt, ns, ymax, betaF, c)
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
  
  results = zeros(totalNS, 2);
  
  istart = 1;
  % three for block stacked
  for rt0=1:nt
      y0 = ylist(rt0);
      t0 = tlist(rt0);
      sechct = sech(c*t0);
      % dydt = a + b*c*sechct*sechct; % use this
      for rt1=1:nt
        % y1= ylist(rt1); % y1 not used. Only for debug usage.
        sigma01 = sigmatable(rt0, rt1); % also multiply this in entries
        ns0 = NStable(rt0, rt1);
        % Attention. Due to the optimization, s0 has to loop 
        % from small to large.
        % index of nonzero M rows: istart~istart+ns0-1
        for rs0=1:ns0
          s0 = slist(rs0); % use this
          s0elm = sigma01+s0;
          if rs0 == ns0
            % analytical
            s0elm = (1+s0elm*betaF)/betaF;  % <x> outside the regime
          end
          
          results(istart,:) = [y0, s0elm];
        % three for block ended
          istart = istart + 1;
        end
      end
  end
end