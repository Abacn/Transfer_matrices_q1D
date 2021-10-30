function Y = NNN3DMatrix(X, nr, nt, ns, dmax, betaF, option)

% The particle diameter is set to be 1.0.
% This function returns Y = A * X, where X is the input vector and A is the
% transfer matrix of dimension (nr^2 * nt * ns) * (nr^2 * nt * ns).
% dmax is the maximum d that the particle can reach, dmax = D-d).
% nr, nt and ns are the number of grids for rh, theta and x.
% betaF is the pressure multiplied by the inverse temperature.
% The advantage of this function is that one does not have to store the
% huge matrix A in the calculation of A * X because A(i, j) can be
% calculated by i and j.

% option: size or calc (default)
% size: return the required vector size

% *** Remind that drh*dt/2*exp(-betaF*xmin) is not included in the matrix ***
% *** dy = ymax/(ny-1); ds see code. ***

global gcount;

if nargin<7
  option = 'calc';
end
% betaP = betaF  / (pi*((dmax+1)/2)^2);
rhmax = (dmax/2)^2;
deltarh = rhmax/nr;
rhlist = linspace(rhmax-deltarh, 0, nr) + deltarh/2; % midpoint rule
deltat = 2*pi/nt;
tlist = linspace(0, 2*pi-deltat, nt);
if dmax < 1.0
    xmin = sqrt(1.0 - dmax^2);
else
    xmin = 0.0;
end
smaxmax = 1.0 - xmin;
deltas = smaxmax / (ns-1);
slist = linspace(0, smaxmax, ns);

btable = zeros(nt, nr, nr);      % table of Boltzmann weights
sigmatable = zeros(nt, nr, nr);  % table of sigma

for rp=1:nr
  for rq=rp:nr
    rh0 = rhlist(rp);
    rh1 = rhlist(rq);
    sigmatmp = 1-rh0-rh1+2*sqrt(rh0*rh1)*cos(tlist);
    sigmatmp(sigmatmp < 0) = 0;
    sigmalist = sqrt(sigmatmp);
    blist = exp(-betaF*(sigmalist-xmin));
    btable(:, rp, rq) = blist;
    sigmatable(:, rp, rq) = sigmalist;
    if rp ~= rq
      btable(:, rq, rp) = blist;
      sigmatable(:, rq, rp) = sigmalist;
    end
  end
end
NStable = ceil((1.0-sigmatable)/deltas)+1;
NStable(NStable > ns) = ns;
totalNS = sum(NStable(:));
lenX = length(X);
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
for rrh1=1:nr
  joffset = 0;
  for rrh2=1:nr
    for rt1=1:nt
      kstart = jstart+joffset;
      ns1=NStable(rt1, rrh1, rrh2);
      nowk = kstart+ns1-1;
      tmpsum = 0.0;
      for irs0=nowk:-1:kstart
        tmpsum = tmpsum + X(irs0);
        Stable(irs0) = tmpsum;
      end
      joffset = joffset + ns1;
    end
  end
  jstart = jstart + joffset;
end

istart = 1;
for rrh0=1:nr
  jstart = 1;
  for rrh1=1:nr
    NS2table = NStable(:,:,rrh1);
    for rt0=1:nt
      ns0=NStable(rt0, rrh0, rrh1);
      sigma01 = sigmatable(rt0, rrh0, rrh1);
      NNNendtable = NS2table;      % index where the NNN interaction ends
      allflag = false;             % the whole row ends NNN interaction
      for rs0=1:ns0
        s0 = slist(rs0);
        pre_bzm = btable(rt0, rrh0, rrh1)*exp(-betaF*s0);
        if rs0==ns0
          pre_bzm = pre_bzm/(betaF*deltas);
        end
        if ~allflag
          allflag = true;
          joffset = 0;
          for rrh2=1:nr
            for rt1=1:nt
              NNNendind = NNNendtable(rt1, rrh2);
              kstart = jstart+joffset;
              if NNNendind ~= 0
                sigma12 = sigmatable(rt1, rrh2, rrh1);
                sigma02 = sigmatable(mod(rt0+rt1-2, nt)+1, rrh2, rrh0);
                s1min = sigma02 - sigma01 - sigma12 - s0;
                if s1min > 0.0 % particle 0 and 2 may overlap
                  soffset = round(s1min/deltas);
                  if soffset >= NNNendind
                    error('Boom! Unexpected overlapped.')
                  end
                  allflag = false;
                  Y(istart) = Y(istart) + Stable(kstart+soffset);
                else
                  Y(istart) = Y(istart) + Stable(kstart);
                  NNNendtable(rt1, rrh2) = 0;
                end
              else
                Y(istart) = Y(istart) + Stable(kstart);
                %Y(istart) = Y(istart) + sum(X(kstart:kstart+NNNendind-1));
              end
              joffset = joffset + NS2table(rt1, rrh2);
            end
          end
          if allflag
            NNNsto = Y(istart);
          end
        else
          Y(istart) = NNNsto;
        end % if all flag
        Y(istart) = Y(istart)*pre_bzm;
        istart = istart + 1;
      end % rs0
    end % rrt0
    jstart = jstart + joffset;
  end % rrh1
end % rrh0

gcount = gcount + 1;
fprintf('calling count: %d\n', gcount);

end
% Nx =  NNN3DMatrix(0, 25, 25, 25, 0.5, log(2), 'size'); X = ones(Nx,1);Y = NNN3DMatrix(X, 25, 25, 25, 0.5, log(2));
% Nx =  NNN3DMatrix(0, 50, 50, 50, 1.0, 5.0, 'size'); X = ones(Nx,1);Y = NNN3DMatrix(X, 50, 50, 50, 1.0, 5.0);
% lambda1 = eigs(@(X) NNN2DMatrix(X, 100, 100, 0.7, 7.0, 3), Nx, 3, 'LM')