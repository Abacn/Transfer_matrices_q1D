% Calculate xi from g(|i-j|)

function xis = redoxis()
    %folder = fileparts(which(mfilename));
    %addpath(strcat(folder, '/../MATLAB'));
    gs = dlmread('gs.dat','',0,0);
    gs = gs(2:end, :);
    nxis = size(gs, 2);
    uplimit = size(gs, 1);
    allfitydata = log(abs(gs));
    for rq=1:nxis
        for rightc=1:uplimit
            if allfitydata(rightc,rq) < log(1e-12) % very faint correlation
                break;
            end
        end
        leftc = max(2, rightc-4);
        if rightc-leftc<4   % very short correlation
            rightc = leftc+4;
        end

        fitxdata = (leftc:rightc).'-1;
        fitydata = allfitydata(leftc:rightc, rq);
        [b, bint] = regress(fitxdata, [fitydata, ones(numel(fitydata),1)]);
        berr = abs((bint(1,2)-bint(1,1))/b(1))*0.5;
        if berr > 0.05
            leftc = min(4, leftc); % need to fix. Try increase fitting range
            fitxdata = (leftc:rightc).'-1;
            fitydata = allfitydata(leftc:rightc, rq);
            [b, bint] = regress(fitxdata, [fitydata, ones(numel(fitydata),1)]);
        end
        xis(rq,:) = [-b(1), -bint(1,2), -bint(1,1)];
    end
    dlmwrite('xis.dat', xis, 'delimiter', '\t', 'Precision', 10);
end