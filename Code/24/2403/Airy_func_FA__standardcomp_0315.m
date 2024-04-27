clc
clear

Nspace = 140;%real N = N+1
[D,x] = cheb(Nspace);

N = Nspace;

D2 = D^2;
D2 = D2(2:N,2:N);

% D2 = D2(1:N+1,1:N+1);

[V,lam] = eig(D2,diag(x(2:N)));
% [V,lam] = eig(D2,diag(x(1:N+1)));

lam = diag(lam);

ii = find(lam>0);

V = V(:,ii);
lam = lam(ii);

% plot(real(lam),'o')
% hold on

[LAM,idx] = sort(lam);
%%
plot(real(LAM),'o')
hold on
%%
Vnew = [zeros(1,length(LAM));V(:,idx);zeros(1,length(LAM))];
Vnew = Vnew/Vnew(N/2+1)*airy(0); 
% figure()
plot(x,Vnew(:,5),'-o')
hold on

%%

% standard dif method

Nreal = N+1;
LHS = zeros(Nreal);
RHS = zeros(Nreal);

Xs = linspace(-1,1,Nreal);
h = 2/N;

K3 = 1/h^2;
K2 = -2/h^2;
K1 = 1/h^2;

for i = 1:Nreal
    for j = 1:Nreal
        if i == 1
            LHS(1,1) = 1;
        elseif i == Nreal
            LHS(Nreal,Nreal) = 1;
        else
            LHS(i,i) = K2;
            LHS(i,i-1) = K3;
            LHS(i,i+1) = K1;
        end
    end
end
RHS(2:Nreal-1,2:Nreal-1) = diag(Xs(1,2:Nreal-1));

[VV,DD] = eigs(LHS,RHS,141,'smallestabs');
Omegalist = diag(DD);

RR = find(Omegalist>0);

VV = VV(:,RR);
Omegalist = Omegalist(RR);

% plot(Omegalist,'o')
% hold on
plot(Xs,-VV(:,5),'-o')
hold on

%%
% ideal airy func
Xideal = -1:0.01:1;
AIRY = airy(7.94413359.*Xideal);

plot(Xideal,AIRY)
hold on

%%
% CHEB  compute D = differentiation matrix, x = Chebyshev grid

  function [D,x] = cheb(N)
  if N==0, D=0; x=1; return, end
  x = cos(pi*(0:N)/N)'; 
  % x = -cos(pi * (0:N)/N)';
  c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
  X = repmat(x,1,N+1);
  dX = X-X';                  
  D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
  D  = D - diag(sum(D'));                 % diagonal entries
  end