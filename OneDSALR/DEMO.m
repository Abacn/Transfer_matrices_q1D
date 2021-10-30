% Transfer matrix method for 1D continuous space SALR model
% Author:  @Abacn

%% Demonstration of the usage of the functions
close all;

% force field coefficients
coeffs = [1, 2.5, 4, 1, 1];

%% Draw a plot of density equation of state of T=0.2
fprintf('Calculate equation of state...\n');
% Temperature
T = 0.2;
% Pressure
ps = [logspace(-6,-1, 20)];
betaPs = ps/T;
% store density
rhos = zeros(1,length(ps));

for rp=1:length(ps)
    p = ps(rp);
   fprintf('P=%.2e\n',p);
   % divide=150 for faster demonstration
   rhos(rp) = findrho(p, 1/T, coeffs, 150);
end
hs = (betaPs-rhos) ./ (rhos.^2);
fprintf('Plot rho over beta*P\n\n');
figure('Position', [200, 300, 800, 600]);
subplot(2,2,1);
loglog(rhos, hs);
xlim([1e-5, 1e-1]);
xlabel('$\rho$', 'Interpreter', 'Latex');
ylabel('$h(\rho)$', 'Interpreter', 'Latex');
title('Equation of state')

%% Draw cluster distribution function K(n) of T=0.2, rho=0.01
fprintf('Calculate cluster distribution function...\n');
T = 0.2;
rho = 0.001;
p = findp(rho, T, coeffs);
ks = cdf(p, 1/T, coeffs);

fprintf('Plot cluster distribution function\n\n');
subplot(2,2,2);
plot(1:length(ks), ks, '-xk');
xlabel('$n$', 'Interpreter', 'Latex');
ylabel('$K(n)$', 'Interpreter', 'Latex');
xlim([1 10]);
ylim([0 1]);
title('Cluster distribution function (CDF)');

%% Draw gap distribution function P_gap(s)
fprintf('Calculate gap distribution function...\n');
[rlist, pros] = gapdf(p, 1/T, coeffs);

fprintf('Plot gap distribution function\n\n');
subplot(2,2,3);
% Remove the last element which is the ideal gas regime s.t. s>\kappa \sigma
plot(rlist(1:end-1), pros(1:end-1));
xlabel('$s$', 'Interpreter', 'Latex');
ylabel('$P_\mathrm{gap}(s)$', 'Interpreter', 'Latex');
title('Gap distribution function (GDF)');

%% Draw correlation lengths for systems with rho=0.01
fprintf('Calculate correlation lengths. This may take for a while...\n');
Ts = 0.05:0.05:0.5;
rho = 0.01;
cors = zeros(1,length(Ts));

for rp=1:length(Ts)
   fprintf('T=%.2f\n', Ts(rp));
   p = findp(rho, Ts(rp), coeffs, 150);
   [~, D] = corlen(p, 1/Ts(rp), coeffs, 150);
   diagD = sort(abs(diag(D)), 'descend');
   cors(rp) = 1/log(diagD(1)/diagD(2));
end
subplot(2,2,4);
plot(Ts, cors);
xlabel('$T$',  'Interpreter', 'Latex');
ylabel('$\xi_s$',  'Interpreter', 'Latex');
title('Correlation length');