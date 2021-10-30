function result = NN2DGr(param, dr, rmax)
    %% getg(r) from sofq
    nargin = 0;
    if nargin==2
        rmax = dr;
        dr = param;
        load param.dat
    elseif nargin==0
        load param.dat
        if nargin<2
            warning('nargin<2');
            dr = 0.05;
            rmax = 180;
        end
    end

    kmax = 2*pi/dr;
    dk = pi/rmax;
    % debug
    rho = NN2DDensity(param(1), param(2), param(3), param(4));
    kresult = NN2DSq(param, dk, kmax);
    kresult(1,2) = 0;
    
    invft = ifft(kresult(:,2) - 1);
    gx = 1 + real(invft)*2*kmax/(2*pi*rho);
    % fix zero frequency shifts
    %invft = invft - (invft(floor(end/2)) - 1);
    % length is dk*kmax
    xs = (0:(length(gx)-1)).'*2*pi/(kmax);
    if 0==nargout
        close all;
        plot(xs(1:floor(end/2)), gx(1:floor(end/2))-1);
        set(gca, 'xscale', 'log', 'yscale', 'log');
        xlim([0.5, rmax]); % xs(end)/2
        ylim([1e-3,3e1]);
        xlabel('$x$', 'interpreter', 'latex');
        ylabel('$g(x)-1$', 'interpreter', 'latex');
        set(gca, 'fontname', 'times new roman', 'fontsize', 18);
        
        hold on;
        fitxs = 10.^linspace(0.01,log10(200));
        ks = (fitxs ./ xs0 + 1 / phat)/(1+1/phat);
        fitgs = cnst*(1+phat)./sqrt(2*pi*(ks-1)) - 1;
        plot(fitxs, fitgs, '-.r');
    end
end