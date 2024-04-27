clc
clear
% 用来算平凡解

%系数
% A=u_a^*
% B=\frac{1}{\rho\phi}
% C=\alpha
% D=k_a
% E=k_{\theta}
% F=k_b
% H=\frac{1}{h_2}
% X=\lambda - here is F
A=1;
B=1;
C=pi/6;

F=0;

H=1;

%网格
L = 0.1;
x = linspace(-L,L,100);
y = linspace(-100*L,100*L,100);

[X,Y]=meshgrid(x,y);

%解的函数关系
UP=A.*(-2.*X.^5 - 4.*X.^3.*Y.^2.*H.^2 - 4.*X.^3.*F.^2 - X.^3.*H.^2.*cos(2.*C) - X.^3.*H.^2 - 2.*X.^2.*F.*H.^2.*sin(2.*C) - 2.*X.*Y.^4.*H.^4 - 4.*X.*Y.^2.*F.^2.*H.^2 + X.*Y.^2.*H.^4.*cos(2.*C) + X.*Y.^2.*H.^4 - 2.*X.*F.^4 + X.*F.^2.*H.^2.*cos(2.*C) - X.*F.^2.*H.^2 + Y.^2.*F.*H.^4.*sin(2.*C))./(2.*(H.^2.*(X.*cos(C) + F.*sin(C)).^2 + (X.^2 + Y.^2.*H.^2 + F.^2).^2));
%
subplot(2,2,1)
surf(X,Y,UP);

zlim([-inf,inf])

legend('Im[\lambda] at k_b=0 of real-state','Location','northwest')

%colorbar;
xlabel('k_a')
ylabel('k_{\theta}')
zlabel('Im[\lambda]')
title('k_b=0')

%
F=1;

%解的函数关系
UP=A.*(-2.*X.^5 - 4.*X.^3.*Y.^2.*H.^2 - 4.*X.^3.*F.^2 - X.^3.*H.^2.*cos(2.*C) - X.^3.*H.^2 - 2.*X.^2.*F.*H.^2.*sin(2.*C) - 2.*X.*Y.^4.*H.^4 - 4.*X.*Y.^2.*F.^2.*H.^2 + X.*Y.^2.*H.^4.*cos(2.*C) + X.*Y.^2.*H.^4 - 2.*X.*F.^4 + X.*F.^2.*H.^2.*cos(2.*C) - X.*F.^2.*H.^2 + Y.^2.*F.*H.^4.*sin(2.*C))./(2.*(H.^2.*(X.*cos(C) + F.*sin(C)).^2 + (X.^2 + Y.^2.*H.^2 + F.^2).^2));
%
subplot(2,2,2)
surf(X,Y,UP);

zlim([-inf,inf])

legend('Im[\lambda] at k_b=1 of real-state','Location','northwest')

%colorbar;
xlabel('k_a')
ylabel('k_{\theta}')
zlabel('Im[\lambda]')
title('k_b=1')

%
F=5;

%解的函数关系
UP=A.*(-2.*X.^5 - 4.*X.^3.*Y.^2.*H.^2 - 4.*X.^3.*F.^2 - X.^3.*H.^2.*cos(2.*C) - X.^3.*H.^2 - 2.*X.^2.*F.*H.^2.*sin(2.*C) - 2.*X.*Y.^4.*H.^4 - 4.*X.*Y.^2.*F.^2.*H.^2 + X.*Y.^2.*H.^4.*cos(2.*C) + X.*Y.^2.*H.^4 - 2.*X.*F.^4 + X.*F.^2.*H.^2.*cos(2.*C) - X.*F.^2.*H.^2 + Y.^2.*F.*H.^4.*sin(2.*C))./(2.*(H.^2.*(X.*cos(C) + F.*sin(C)).^2 + (X.^2 + Y.^2.*H.^2 + F.^2).^2));
%
subplot(2,2,3)
surf(X,Y,UP);

zlim([-inf,inf])

legend('Im[\lambda] at k_b=5 of real-state','Location','northwest')

%colorbar;
xlabel('k_a')
ylabel('k_{\theta}')
zlabel('Im[\lambda]')
title('k_b=5')

%
F=10;

%解的函数关系
UP=A.*(-2.*X.^5 - 4.*X.^3.*Y.^2.*H.^2 - 4.*X.^3.*F.^2 - X.^3.*H.^2.*cos(2.*C) - X.^3.*H.^2 - 2.*X.^2.*F.*H.^2.*sin(2.*C) - 2.*X.*Y.^4.*H.^4 - 4.*X.*Y.^2.*F.^2.*H.^2 + X.*Y.^2.*H.^4.*cos(2.*C) + X.*Y.^2.*H.^4 - 2.*X.*F.^4 + X.*F.^2.*H.^2.*cos(2.*C) - X.*F.^2.*H.^2 + Y.^2.*F.*H.^4.*sin(2.*C))./(2.*(H.^2.*(X.*cos(C) + F.*sin(C)).^2 + (X.^2 + Y.^2.*H.^2 + F.^2).^2));
%
subplot(2,2,4)
surf(X,Y,UP);

zlim([-inf,inf])

legend('Im[\lambda] at k_b=10 of real-state','Location','northwest')

%colorbar;
xlabel('k_a')
ylabel('k_{\theta}')
zlabel('Im[\lambda]')
title('k_b=10')





%title('characteristic plane Im[\lambda]-(k_a,k_{\theta}) at special k_b of real-state','FontSize',30)





