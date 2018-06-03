function [xi, V, D] = NN2DCalc(param)
    if nargin<1
        % test
        load param.dat;
    end
    ny1 = param(1);
    ymax = param(2);
    betaP = param(3);
    tanhc = param(4);

    [V, D] = NN2DEigs(ny1, ymax, betaP, tanhc);

    xi(1) = 1.0 / log(abs(D(1, 1)) / abs(D(2, 2)));
    xi(2) = 1.0 / log(abs(D(1, 1)) / abs(D(3, 3)));
end

% [xi, V, D]=NN2DCalc();
