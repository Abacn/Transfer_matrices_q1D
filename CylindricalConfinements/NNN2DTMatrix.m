function Y = NNN2DTMatrix(X, nt, ns, ymax, betaF, c, option)

% Same as NNN2DMatrix, but it returns Y = A.T * X, or Y.T = x.T*A
% where X is the input vector and A is the transfer matrix
% basically this function is used to get the left eigenvector

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
NStable = ceil((1.0-sigmatable)/deltas)+1; % 1e-6: To avoid numerical err
NStable(NStable > ns) = ns;
% total vector length
totalNS = sum(NStable(:));
if lenX ~= totalNS
    if strcmp(option, 'calc')
        error(sprintf('Vector X size should be %d', totalNS));
    else
        Y = totalNS;
      return 
    end
end
Y = zeros(lenX, 1);

% indices of rearranging rows
raind = rearrangerows(NStable);
raX = X(raind);
Stable = zeros(lenX, 1);
% build sumtable
iend = 1;
for rt1=1:nt
  for rt0=1:nt
    t0 = tlist(rt0);
    sechct = sech(c*t0);
    dydt = a + b*c*sechct*sechct;
    ns0=NStable(rt0, rt1);
    tmpsum = 0.0;
    sigma0 = dydt*btable(rt0, rt1);
    % sigma01 = sigmatable(rt0, rt1);
    iend = iend+ns0;
    for irs0=1:ns0
      nowk = iend-irs0;
      % pre_bzm = dydt*exp(-betaF*(+sigma01-xmin+slist(ns0+1-irs0)));
      pre_bzm = sigma0*exp(-betaF*slist(ns0+1-irs0));
      if irs0==1
        pre_bzm = pre_bzm/(betaF*deltas);
      end
      tmpsum = tmpsum + raX(nowk)*pre_bzm;
      Stable(nowk) = tmpsum;
    end
  end
end

istart = 1;
jstart = 1;
% three for block stacked
for rt1=1:nt
  NS2table = NStable(:,rt1);
  for rt2=1:nt
    ns1 = NS2table(rt2);
    sigma12 = sigmatable(rt1, rt2); % also multiply this in entries
    NNNendtable = NStable(:,rt1);   % s-wise flag
    allflag = false;             % the whole row ends NNN interaction
    % index of nonzero M rows: istart~istart+ns0-1
    for rs1=1:ns1
      s1 = slist(rs1); % use this
      if ~allflag
        allflag = true;
        joffset = 0;
        for rt0=1:nt
          NNNendind = NNNendtable(rt0);
          kstart = jstart+joffset;
          if NNNendind ~= 0
            sigma01 = sigmatable(rt0, rt1);
            sigma02 = sigmatable(rt0, rt2);
            s0min = sigma02 - sigma01 - sigma12 - s1;
            if s0min > 0.0 % particle 0 and 2 may overlap
              soffset = round(s0min/deltas);
              if soffset >= NNNendind
                  error('Boom! unexpected contact happened.')
              end
              allflag = false;
              Y(istart) = Y(istart) + Stable(kstart+soffset);
            else
              Y(istart) = Y(istart) + Stable(kstart);
              NNNendtable(rt0) = 0;
            end
          else
            Y(istart) = Y(istart) + Stable(kstart);
          end
          joffset = joffset + NS2table(rt0);
        end
        if allflag
          NNNsto = Y(istart);
        end
      else
        Y(istart) = NNNsto;
      end % if all flag
      istart = istart + 1;
    end % rs0
  end % rt2
  jstart = jstart + joffset;
end % rt1

gcount = gcount + 1;
fprintf('calling count: %d\n', gcount);
end

% return a indices list of mapping (t0, t1, s0) to be (t1, t0, s0)
function indices = rearrangerows(NStable)
  nt = size(NStable, 1);
  indtable = zeros(nt, nt);
  lastind = 1;
  for rt0=1:nt
    for rt1=1:nt
      indtable(rt0, rt1) = lastind;
      lastind = lastind + NStable(rt0, rt1);
    end
  end
  totalNS = sum(NStable(:));
  indices=zeros(totalNS, 1);
  lastind = 1;
  for rt1=1:nt
    for rt0=1:nt
      ns0 = NStable(rt0, rt1);
      istart = indtable(rt0, rt1);
      indices(lastind:lastind+ns0-1) = (istart:istart+ns0-1).';
      lastind = lastind+ns0;
    end
  end
end

% Nx =  NNN2DMatrix(0, 100, 100, 0.7, 7.0, 3, 'size'); X = ones(Nx,1);Y = NNN2DMatrix(X, 100, 100, 0.7, 7.0, 3);
% lambda1 = eigs(@(X) NNN2DMatrix(X, 100, 100, 0.7, 7.0, 3), Nx, 3, 'LM')