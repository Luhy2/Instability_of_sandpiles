% p15.m - solve eigenvalue BVP u_xx = lambda*u, u(-1)=u(1)=0
clc
clear
%   N = 36;
N=1000;
  [D,x] = cheb(N);
  D2 = D^2;
  D2 = D2(2:N,2:N);
  
  [V,Lam] = eig(D2);
  lam = diag(Lam);
  [foo,ii] = sort(-lam);          % sort eigenvalues and -vectors
  
  lam = lam(ii); V = V(:,ii);
  TRY=280;
  for j = TRY:TRY                  % plot 6 eigenvectors
     u = [0;V(:,j);0]; %subplot(7,1,j)
%     plot(x,u,'.','markersize',12), grid on
plot(x,u,'-o'), grid on
hold on
%     xx = -1:.01:1; uu = polyval(polyfit(x,u,N),xx);
%     line(xx,uu), axis off
    
  end
  xc=linspace(-1,1,10*N);
 COMP=-0.046*sin(xc.*sqrt(-lam(TRY,1)));
 COMP2=1*cos(xc.*sqrt(-lam(TRY,1)));
 plot(xc,-COMP,'LineWidth', 1)
 hold on

 legend
  %%
  Length=length(lam);
  LST=linspace(1,Length,Length);
  
  plot(LST,-lam,'o')
  hold on
  
sqrtk=((LST+0.5).*pi)./2;
doubk=sqrtk.^2;
plot(LST,doubk,'o')
hold on

sqrtk2=((LST).*pi)./2;
doubk2=sqrtk2.^2;
plot(LST,doubk2,'o')
hold on

  
  %%
  
  
  
% CHEB  compute D = differentiation matrix, x = Chebyshev grid

%   function [D,x] = cheb(N)
%   if N==0, D=0; x=1; return, end
%   x = cos(pi*(0:N)/N)'; 
%   c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
%   X = repmat(x,1,N+1);
%   dX = X-X';                  
%   D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
%   D  = D - diag(sum(D'));                 % diagonal entries
%   end
  
  
  