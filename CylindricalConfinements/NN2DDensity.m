function rho = NN2DDensity(nt, ymax, betaF, c, V, D)
    % Called by script, output to the working dir
    % two mode
    % (1) parameter input by argin: NN2DDensity(nt, ns, betaF, tanhc, (optional: V, D))
    % (2) parameter input by file param.dat: NN2DDensity((optional: V, D))
    % if V, D is assigned, used it to get qM and eigM, save half of time.
    % flag for eigenvalue and eigenvector input
    egflg = false;
    saveflag = false;
    
    if nargin==2
      qM = nt(:,1);
      eigM = ymax(1,1)/prefix;
      egflg = true;
    elseif nargin==6
      xmin = sqrt(1.0 - ymax*ymax);
      deltat = 2.0 / nt;
      prefix = exp(-betaF*xmin)*deltat;
      qM = V(:,1);
      eigM = D(1,1)/prefix;
      egflg = true;
    end
    if nargin<4 
      load param.dat;
      nt = param(1);
      ymax = param(2);
      betaF = param(3);
      c = param(4);
    end
    if betaF<0
        betaF = abs(betaF);
    end
    deltat = 2.0 / (nt-1);
    xmin = sqrt(1.0 - ymax*ymax);
    prefx = exp(-betaF*xmin)*deltat;
    reloadflag =false;
    stofname = 'NN2DDensity_redo';
    
    if saveflag && exist(strcat(stofname,'.mat'), 'file')
        load(stofname);
        if param_sv(1)==nt && param_sv(2)==ymax && param_sv(3)==betaF && param_sv(4)==c
            reloadflag = true;
        end 
    end
    
    
    if ~(egflg || reloadflag)
      fprintf('Calculating q\n');
      [qM, eigM] = NN2DEigs(nt, ymax, betaF, c, 1, false);
      % test: finite consistency, differentiation
      %deltabetaF = 0.0001;
      %[qMb, eigMb] = NN2DEigs(nt, ymax, betaF*(1+deltabetaF), c, 1, true);
      %rho = betaF/ ( 1.0 + log(eigM/eigMb)./deltabetaF);
      if saveflag
        param_sv = [nt, ymax, betaF, c];
        save(stofname, 'qM', 'eigM', 'param_sv');
      end
    end
    
    rho = 1.0/(1.0/betaF + qM.' * NN2DDifmat(nt,ymax, betaF, c) * qM / eigM(1,1));
    %dlmwrite('rho.dat', rho, 'delimiter', '\t', 'Precision', 10);

    param_sv = [nt, ymax, betaF, c];
    save('NN2Drho', 'qM', 'eigM', 'param_sv', 'rho');
    % test
    %{
    b = ymax / 2.0;
	a = ymax / 2.0 - b * tanh(c);
    dt = 2.0/nt;
	tlist = linspace(-1+deltat/2, 1-deltat/2, nt);
	ylist = tlist*a + tanh(tlist*c)*b;
    sechct = sech(c*tlist);
    dydt = a + b*c*sechct.^2; % dy/dt
    normterm = dydt * (qM.^2);
    %plot(tlist, (qM.^2) );
    plot(ylist, (qM.^2)/normterm );
    %}
end
% rho = NNN2DDensity(100, 100, 0.5, 5, 3);


% return d K(\Delta x) /d (\beta F), results.*K gives the tangent matrix
function results = NN2DDifmat(nt, ymax, betaF, c)
	b = ymax / 2.0;
	a = ymax / 2.0 - b * tanh(c);

	deltat = 2.0 / nt;
    xmin = sqrt(1.0 - ymax*ymax);
	tlist = linspace(-1+deltat/2, 1-deltat/2, nt);
	ylist = tlist*a + tanh(tlist*c)*b;
	tmp1 = repmat(ylist,nt,1); tmp2 = tmp1-tmp1.'; tmp3 = tmp2.*tmp2; tmp3(tmp3>1.0)=1.0;
	sigmatable = sqrt(1.0-tmp3);
	sechct = sech(c*tlist);
    sqdydt = sqrt(a + b*c*sechct.^2); % sqrt(dy/dt)
    tmp4 = sqdydt.' * sqdydt;
	results = exp(-betaF*(sigmatable-xmin)) .* sigmatable .* tmp4;
end
