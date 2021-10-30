function Y = NNN2DMatrix(X, nt, ns, ymax, betaF, c, option)

% The particle diameter is set to be 1.0.
% This function returns Y = A * X, where X is the input vector and A is the
% transfer matrix of dimension (nt * nt * ns) * (nt * nt * ns).
% ymax is the maximum y that the particle can reach, ymax = H-d).
% nt and ns are the number of grids for y1, y2 and x.
% betaF is the pressure multiplied by the inverse temperature.
% The advantage of this function is that one does not have to store the
% huge matrix A in the calculation of A * X because A(i, j) can be
% calculated by i and j.

% option: size or calc (default)
% size: return the required vector size

% *** Remind that dt*ds*exp(-betaF*xmin) is not included in the matrix ***
% *** dy = ymax/(ny-1); ds see code. ***

global gcount;

if nargin<7
  option = 'calc';
end
nt21 = nt*nt;
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
lenX = length(X);

% two end all considered. May test if midpoint rule results sth. different.
tlist = linspace(-1, 1, nt);
slist = linspace(0, smaxmax, ns); % scaled s. The first element denotes smax, and it is ideal
ylist = tlist*a + tanh(tlist*c)*b;
tmp1 = repmat(ylist,nt,1); tmp2 = tmp1-tmp1.'; tmp3 = tmp2.*tmp2; tmp3(tmp3>1.0)=1.0;
sigmatable = sqrt(1.0-tmp3);
btable = exp(-betaF*(sigmatable-xmin));
% the number of s we have
NStable = ceil((1.0-sigmatable)/deltas)+1;
NStable(NStable > ns) = ns;
% total vector length
totalNS = sum(NStable(:));
if lenX ~= totalNS
    if strcmp(option, 'calc')
        error('Vector X size should be %d', totalNS);
    else
        Y = totalNS;
      return 
    end
end

Y = zeros(lenX, 1);
Stable = zeros(lenX, 1);
% build sumtable
jstart = 1;
for rt1=1:nt
  joffset = 0;
  for rt2=1:nt
    kstart = jstart+joffset;
    ns1=NStable(rt1, rt2);
    nowk = kstart+ns1-1;
    tmpsum = 0.0;
    for irs0=nowk:-1:kstart
      tmpsum = tmpsum + X(irs0);
      Stable(irs0) = tmpsum;
    end
    joffset = joffset + ns1;
  end
  jstart = jstart + joffset;
end

istart = 1;
% three for block stacked
for rt0=1:nt
  % y0 = ylist(rt0); y0 not used. Only for debug usage.
  t0 = tlist(rt0);
  sechct = sech(c*t0);
  dydt = a + b*c*sechct*sechct; % use this
  jstart = 1;
  for rt1=1:nt
    % y1 = ylist(rt1);
    ns0 = NStable(rt0, rt1);
    % NS2table = NStable(:,rt1);
    sigma0 = dydt*btable(rt0, rt1);
    sigma01 = sigmatable(rt0, rt1); % also multiply this in entries
    NNNendtable = NStable(:,rt1);      % index where the NNN interaction ends
    allflag = false;             % the whole row ends NNN interaction
    % Attention. Due to the optimization, s0 has to loop 
    % from small to large.
    % index of nonzero M rows: istart~istart+ns0-1
    for rs0=1:ns0
      s0 = slist(rs0); % use this
      pre_bzm = sigma0*exp(-betaF*s0);
      % pre_bzm = exp(-betaF*(s0+sigma01-xmin))*dydt;
      if rs0 == ns0
        % analytical
        pre_bzm = pre_bzm/(betaF*deltas);
      end
      if ~allflag
        allflag = true;
        joffset = 0;
        for rt2=1:nt
          ns1 =  NStable(rt1, rt2);
          NNNendind = NNNendtable(rt2);
          % index of start column: jstart~jstart+ns
          kstart = jstart+joffset;
          if NNNendind ~= 0
            sigma12 = sigmatable(rt1, rt2);
            sigma02 = sigmatable(rt0, rt2);
            s1min = sigma02 - sigma01 - sigma12 - s0;
            if s1min > 0.0 % particle 0 and 2 may overlap
              soffset = round(s1min/deltas);
              allflag = false;
              if soffset >= ns1
                  error('Boom! unexpected contact happened.')
              end
              Y(istart) = Y(istart) + Stable(kstart+soffset);
            else
              Y(istart) = Y(istart) + Stable(kstart);
              NNNendtable(rt2) = 0;
            end
          else
            Y(istart) = Y(istart) + Stable(kstart);
          end
          joffset = joffset + NStable(rt1, rt2);
        end % rt2
        if allflag
          NNNsto = Y(istart);
        end
      else
        Y(istart) = NNNsto;
      end % if all flag
      Y(istart) = Y(istart)*pre_bzm;
      istart = istart + 1;
    end % rs0
    jstart = jstart + joffset;
  end % rt1
end % rt0
      
gcount = gcount + 1;
fprintf('calling count: %d\n', gcount);
end

% Nx =  NNN2DMatrix(0, 100, 100, 0.7, 7.0, 3, 'size'); X = ones(Nx,1);Y = NNN2DMatrix(X, 100, 100, 0.7, 7.0, 3);
% lambda1 = eigs(@(X) NNN2DMatrix(X, 100, 100, 0.7, 7.0, 3), Nx, 3, 'LM')