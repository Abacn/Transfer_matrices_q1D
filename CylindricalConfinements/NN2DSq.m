function result = NN2DSq(param, dk, kmax)
    % calculating structural factor
    % param: parameters
    % dk: discretization gap
    % kmax: max wavenumber
    % Reference: J. F. Robinson, M. J. Godfrey, and M. A. Moore Phys. Rev. E 93, 032101
    if nargin==2
        kmax = dk;
        dk = param;
        load param.dat
    elseif nargin==0
        load param.dat
        if nargin<2
            warning('nargin<2');
            dk = 0.05;
            kmax = 30;
        end
    end
    nt = param(1);
    ymax = param(2);
    betaF = param(3);
    c = param(4);
    
    fprintf('Calculating q\n');
    [qM, eigM] = NN2DEigs(nt, ymax, betaF, c, 1, true);
    qM = abs(qM);
    ks = (0:dk:kmax).';
    sq = repmat(ks, 1);
    % solve linear equations
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
    
    for rp=1:length(ks)
        pf = betaF/(eigM*(betaF - 1i*ks(rp)));
        Smat = exp((1i*ks(rp)-betaF)*sigmatable) .* tmp4 * pf * deltat;
        phis = linsolve(eye(nt)-Smat,qM);
        sq(rp) = 1 + 2*real( qM.'*Smat*phis );
    end
    
    if nargout==0
        close all;
        sq(1) = nan;
        plot(ks, sq);
        xlabel('$k$', 'interpreter', 'latex');
        ylabel('$S(k)$', 'interpreter', 'latex');
        set(gca, 'fontname', 'times new roman', 'fontsize', 18);
    else
        result = [ks, sq];
    end
end

