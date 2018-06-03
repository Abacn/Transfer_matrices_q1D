function Y = NNN3DTMatrix(X, nr, nt, ns, dmax, betaF, option)

% Same as NNN3DMatrix, but it returns Y = A.T * X, or Y.T = x.T*A
% basically this function is used to get the left eigenvector

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
% indices of rearranging rows
raind = rearrangerows(NStable);
raX = X(raind);
Stable = zeros(lenX, 1);
% build sumtable
iend = 1;
for rrh1=1:nr
  for rrh0=1:nr
    for rt0=1:nt
      ns0=NStable(rt0, rrh1, rrh0);
      tmpsum = 0.0;
      sigma0 = btable(rt0, rrh1, rrh0);
      iend = iend+ns0;
      for irs0=1:ns0
        nowk = iend-irs0;
        pre_bzm = sigma0*exp(-betaF*slist(ns0+1-irs0));
        if irs0==1
          pre_bzm = pre_bzm/(betaF*deltas);
        end
        tmpsum = tmpsum + raX(nowk)*pre_bzm;
        Stable(nowk) = tmpsum;
      end
    end
  end
end

istart = 1;
jstart = 1;
for rrh1=1:nr
  NS2table = NStable(:,:,rrh1);
  for rrh2=1:nr
    for rt1=1:nt
      ns1=NStable(rt1, rrh1, rrh2);
      sigma12 = sigmatable(rt1, rrh1, rrh2);
      NNNendtable = NS2table;      % s-wise flag
      allflag = false;             % the whole row ends NNN interaction
      for rs1=1:ns1
        s1 = slist(rs1);
        if ~allflag
          allflag = true;
          joffset = 0;
          for rrh0=1:nr
            for rt0=1:nt
              NNNendind = NNNendtable(rt0, rrh0);
              kstart = jstart+joffset;
              if NNNendind ~= 0
                sigma01 = sigmatable(rt0, rrh0, rrh1);
                sigma02 = sigmatable(mod(rt0+rt1-2, nt)+1, rrh2, rrh0);
                s0min = sigma02 - sigma01 - sigma12 - s1;
                if s0min > 0.0 % particle 0 and 2 may overlap
                  soffset = round(s0min/deltas);
                  if soffset >= NNNendind
                    error('Boom! Unexpected overlapped.')
                  end
                  allflag = false;
                  Y(istart) = Y(istart) + Stable(kstart+soffset);
                else
                  Y(istart) = Y(istart) + Stable(kstart);
                  NNNendtable(rt0, rrh0) = 0;
                end
              else
                Y(istart) = Y(istart) + Stable(kstart);
                %Y(istart) = Y(istart) + sum(X(kstart:kstart+NNNendind-1));
              end
              joffset = joffset + NS2table(rt0, rrh0);
            end
          end
          if allflag
            NNNsto = Y(istart);
          end
        else
          Y(istart) = NNNsto;
        end % if all flag
        istart = istart + 1;
      end % rs1
    end % rrt0
  end % rrh2
  jstart = jstart + joffset;
end % rrh1
gcount = gcount + 1;
fprintf('calling count: %d\n', gcount);

end

% return a indices list of mapping (r0, r1, t0, s0) to be (r1, r0, t0, s0)
function indices = rearrangerows(NStable)
  nt = size(NStable, 1);
  nr = size(NStable, 2);
  indtable = zeros(nt, nr, nr);
  lastind = 1;
  for rrh0=1:nr
    for rrh1=1:nr
      for rt0=1:nt
        indtable(rt0, rrh0, rrh1) = lastind;
        lastind = lastind + NStable(rt0, rrh0, rrh1);
      end
    end
  end
  totalNS = sum(NStable(:));
  indices=zeros(totalNS, 1);
  lastind = 1;
  for rrh1=1:nr
    for rrh0=1:nr
      for rt0=1:nt
        ns0 = NStable(rt0, rrh0, rrh1);
        istart = indtable(rt0, rrh0, rrh1);
        indices(lastind:lastind+ns0-1) = (istart:istart+ns0-1).';
        lastind = lastind+ns0;
      end
    end
  end
end

% Nx =  NNN3DMatrix(0, 25, 25, 25, 0.5, log(2), 'size'); X = ones(Nx,1);Y = NNN3DTMatrix(X, 25, 25, 25, 0.5, log(2));
% Nx =  NNN3DMatrix(0, 20, 20, 20, 0.5, log(2), 'size'); X = zeros(Nx,1);X(2000)=1;Y = NNN3DTMatrix(X, 20, 20, 20, 0.5, log(2));plot(Y);
