% Generating g(r)

close all;
choice = 2;

%% load tmat result
ymax_prop = 0.5;
P_prop = 5;
stryp = strcat(num2str(ymax_prop), '_', num2str(P_prop));

datarootdir = '../data/2D/dists/';
load(strcat('NN2Drho_',stryp,'.mat'));
fnamepref = 'NN2Dngh_';

if 1==choice
    %% generate subsequent points
    
    % parameters
    nt = param_sv(1);
    ymax = param_sv(2);
    betaF = param_sv(3);
    c = param_sv(4);
      
    % get cumulative marginal distribution
    b = ymax / 2.0;
    a = ymax / 2.0 - b * tanh(c);

    deltat = 2.0 /nt;
    tlist = linspace(-1+deltat/2, 1-deltat/2, nt);
    tlistb = [-1, tlist, 1];
    ylistb = tlistb*a + tanh(tlistb*c)*b;
    sechctb = sech(c*tlistb);
    dydtb = a + b*c*sechctb.^2;
    sqdydtb = sqrt(dydtb);
    
    cumuprops = cumsum((qM.').^2);
    tcums = linspace(-1, 1, nt+1);
    cumuprops = [0, cumuprops];
    xmin = sqrt(1.0 - ymax*ymax);

    nowN=0;
    Ndup=400;
    Nshift=0;
    maxx = 2000;
    if ~exist(datarootdir, 'dir')
        mkdir(datarootdir);
    end
    datadir = strcat(datarootdir, stryp);
    if ~exist(datadir, 'dir')
        mkdir(datadir);
    end
    while nowN<Ndup
        disx = 0.0;
        ny = 0;
        % choose t from marginal probability
        chooset = interp1(cumuprops, tcums, rand(), 'pchip');
        choosey = a*chooset + b*tanh(c*chooset);
        fout = fopen(strcat(datadir, '/', fnamepref, num2str(nowN+Nshift),'.dat'), 'w');
        fprintf(fout, '%.6f\t%.6f\n', disx, choosey);
        while disx <= maxx
            x0min = sqrt(1-(b+abs(choosey)).^2);
            distable = choosey-ylistb(2:end-1);
            sigmatable = sqrt(1-distable.^2);
            % condition probability distribution P(y'|y)
            sqdydty = interp1(ylistb, sqdydtb, choosey);
            condprop = exp(-betaF*(sigmatable-x0min)) .* sqdydtb(2:end-1) * sqdydty .* qM.';
            % fix identity value of condprop in extremely high pressure
            cumuconprops = zeros(1, length(condprop)+1);
            rp=1;
            while rp <= length(condprop)
                cumuconprops(rp+1) = cumuconprops(rp) + condprop(rp);
                if cumuconprops(rp+1) == cumuconprops(rp)
                    break
                end
                rp = rp + 1;
            end
            cumuconprops = cumuconprops(1:rp);
            %cumuconprops = [0, cumsum(condprop)];
            targetp = cumuconprops(end)*rand();
            nextt = interp1(cumuconprops, tcums(1:length(cumuconprops)), targetp);
            nexty = a*nextt + b*tanh(c*nextt);
            x0 = sqrt(1-(nexty-choosey).^2);
            x1 = -log(1-rand())/betaF;
            disx = disx + x0 + x1;
            choosey = nexty;
            fprintf(fout, '%.6f\t%.6f\n', disx, choosey);
            %fprintf('%.6f\t%.5f\t%.5f\n', x0, x1, nexty);
        end
        fclose(fout);
        nowN = nowN + 1;
    end
elseif 2==choice
    %% plot g(r)
    plotchoice = 1;
    calcnewchoice = true;
    if 1==plotchoice
        maxedge = 10;
    else
        maxedge = 100;
    end
    if calcnewchoice
        % set after
        halfmax = nan;
        jumpx = nan;
        datadir = strcat(datarootdir, stryp);
        files = dir(datadir);
        nfille = 0;

        if maxedge <= 10
            edges = 0:0.1:maxedge;
        elseif maxedge < 50
            edges = [linspace(0, 10-0.01, 1001), linspace(10, maxedge, 1001)];
        elseif maxedge<1000
            edges = [linspace(0, 10-0.01, 1001), 10:0.01:maxedge];
        else
            edges = [linspace(0, 10-0.01, 1001), 10:0.01:99.995, 100:0.05:maxedge];
        end

        counts = zeros(1, length(edges)-1);

        npart = 0;
        nfile = 0;
        lenA = 0;
        maxlen = length(edges)*1000;
        for rp = 1 : length(files)
            fname = files(rp).name;
            if ~startsWith(fname, fnamepref), continue, end
            fprintf('\r%s', fname);
            data = dlmread(strcat(files(rp).folder, '/', files(rp).name));
            if ~isfinite(halfmax)
                maxx = floor(max(data(:,1)));
                halfmax = maxx - maxedge;
                % take at most 200 particles - if too slow
                jumpx = 1;%max(1, floor(halfmax/maxx*size(data,1) / 200));
            end
            for rq=1:jumpx:length(data)
                if data(rq,1) > halfmax, break, end
                for rr=rq+1:length(data)
                    dis = data(rr,1)-data(rq,1);
                    if dis >= maxedge, break, end
                    lenA = lenA + 1;
                    idx = lowerbound(dis, edges);
                    counts(idx) = counts(idx) + 1;
                end
                npart = npart + 1;
            end
            nfile = nfile + 1;
            %if lenA > maxlen
            %    break;
            %end
        end

        dr = diff(edges);
        centers = ave(edges);
        gr = counts./(dr*rho*npart);

        % save file
        fname = strcat(datarootdir, 'NN2Dgr_',stryp,'_',num2str(maxedge), '.dat' );
        dlmwrite(fname, [centers.', gr.'], 'delimiter', '\t', 'precision', '%.6f');
    else
        fname = strcat(datarootdir, 'NN2Dgr_',stryp,'_',num2str(maxedge), '.dat' );
        data = dlmread(fname);
        centers = data(:,1);
        gr = data(:,2);
    end
    if 1==plotchoice
        plot(centers, gr);
        xlim([0.5,maxedge]);
        ylabel('$g(x)$', 'interpreter', 'latex');
    else
        plot(centers, abs(gr-1));
        set(gca, 'yscale', 'log', 'xscale', 'log');
        xlim([0.5,maxedge]);
        ylim([1e-1,10.^(ceil(log10(max(gr))))]);
        ylabel('$|g(x)-1|$', 'interpreter', 'latex');
    end
    xlabel('$x$', 'interpreter', 'latex');
    
    set(gca, 'fontname', 'times new roman', 'fontsize', 18);
end
