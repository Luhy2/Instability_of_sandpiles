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
UP=A.*H.*(-X.^4 - X.^2.*Y.^2.*H.^2 - 2.*X.^2.*F.^2 - X.^2.*H.^2.*cos(C).^2 - 2.*X.*F.*H.^2.*sin(C).*cos(C) - Y.^2.*F.^2.*H.^2 - F.^4 + F.^2.*H.^2.*cos(C).^2 - F.^2.*H.^2).*cos(C)./(H.^2.*(X.*cos(C) + F.*sin(C)).^2 + (X.^2 + Y.^2.*H.^2 + F.^2).^2);

%
subplot(2,2,1)
surf(X,Y,UP);

zlim([-inf,inf])
legend('Re[\lambda] at k_b=0','Location','northwest');
title('k_b=0')
xlabel('k_a')
ylabel('k_{\theta}')
zlabel('Re[\lambda]')

hold on
%
F=1;

%解的函数关系
UP=A.*H.*(-X.^4 - X.^2.*Y.^2.*H.^2 - 2.*X.^2.*F.^2 - X.^2.*H.^2.*cos(C).^2 - 2.*X.*F.*H.^2.*sin(C).*cos(C) - Y.^2.*F.^2.*H.^2 - F.^4 + F.^2.*H.^2.*cos(C).^2 - F.^2.*H.^2).*cos(C)./(H.^2.*(X.*cos(C) + F.*sin(C)).^2 + (X.^2 + Y.^2.*H.^2 + F.^2).^2);

%
subplot(2,2,2)
surf(X,Y,UP);

zlim([-inf,inf])
legend('Re[\lambda] at k_b=1','Location','northwest');
title('k_b=1')
xlabel('k_a')
ylabel('k_{\theta}')
zlabel('Re[\lambda]')

hold on
%
F=5;

%解的函数关系
UP=A.*H.*(-X.^4 - X.^2.*Y.^2.*H.^2 - 2.*X.^2.*F.^2 - X.^2.*H.^2.*cos(C).^2 - 2.*X.*F.*H.^2.*sin(C).*cos(C) - Y.^2.*F.^2.*H.^2 - F.^4 + F.^2.*H.^2.*cos(C).^2 - F.^2.*H.^2).*cos(C)./(H.^2.*(X.*cos(C) + F.*sin(C)).^2 + (X.^2 + Y.^2.*H.^2 + F.^2).^2);

%
subplot(2,2,3)
surf(X,Y,UP);

zlim([-inf,inf])
legend('Re[\lambda] at k_b=5','Location','northwest');
title('k_b=5')
xlabel('k_a')
ylabel('k_{\theta}')
zlabel('Re[\lambda]')

hold on

%
F=10;

%解的函数关系
UP=A.*H.*(-X.^4 - X.^2.*Y.^2.*H.^2 - 2.*X.^2.*F.^2 - X.^2.*H.^2.*cos(C).^2 - 2.*X.*F.*H.^2.*sin(C).*cos(C) - Y.^2.*F.^2.*H.^2 - F.^4 + F.^2.*H.^2.*cos(C).^2 - F.^2.*H.^2).*cos(C)./(H.^2.*(X.*cos(C) + F.*sin(C)).^2 + (X.^2 + Y.^2.*H.^2 + F.^2).^2);

%
subplot(2,2,4)
surf(X,Y,UP);

zlim([-inf,inf])
legend('Re[\lambda] at k_b=10','Location','northwest');
title('k_b=10')
xlabel('k_a')
ylabel('k_{\theta}')
zlabel('Re[\lambda]')

hold on








%title('characteristic plane Re[\lambda]-(k_a,k_{\theta}) at special k_b of real-state')