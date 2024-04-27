clc
clear
% 用来算平凡解


%系数
A=1;
h=1;
a=30;
kb=1;

%函数关系
L = 1;
x = linspace(-L,L,100);
y = linspace(-100*L,100*L,100);

[X,Y]=meshgrid(x,y);


%F = ((X.^2).*Y./kb.*h.*A-0.5*sin(a).*Y/h./kb+(cos(a))^2-h.*X.*Y+(h^2)*(X.^2))./((h^2).*(X.^2)+(cos(a))^2);
F = -A.*X;

%
surf(X,Y,F);
zlim([-inf,inf])
%axis ([-1 10 -1 10 0 5])
xlabel('k_a','FontSize',30)
ylabel('k_{\theta}','FontSize',30)
zlabel('\lambda','FontSize',30)
title('characteristic plane \lambda-(k_a,k_{\theta}) at special k_b of base-state','FontSize',30)
%%
clc
clear

A=1;
B=1;
C=pi/6;

F=0;

H=1;
%函数关系
L = 10;
x = linspace(-L,L,100);
y = linspace(-100*L,100*L,100);

[X,Y]=meshgrid(x,y);

figure()
JB=- ((F.^2 + Y.^2.*H.^2 + X.^2).*(3.*A.*X.^2.*H.*cos(C) + A.*F.^2.*H.*cos(C) + abs(cos(C)).*real((-(F.^2.*1i + H.*sin(C).*F + X.^2.*1i + H.*cos(C).*X).^2).^(1./2)).*abs(A).*abs(H) + 2.*A.*X.*F.*H.*sin(C)))./(2.*((F.^2 + Y.^2.*H.^2 + X.^2).^2 + (X.*H.*cos(C) + F.*H.*sin(C)).^2)) + ((X.*H.*cos(C) + F.*H.*sin(C)).*(2.*A.*X.^3 + 2.*A.*X.*F.^2 + 2.*A.*X.*Y.^2.*H.^2 - A.*X.*H.^2.*cos(C).^2 + imag((-(F.^2.*1i + H.*sin(C).*F + X.^2.*1i + H.*cos(C).*X).^2).^(1./2)).*abs(cos(C)).*abs(A).*abs(H) - A.*F.*H.^2.*cos(C).*sin(C)))./(2.*((F.^2 + Y.^2.*H.^2 + X.^2).^2 + (X.*H.*cos(C) + F.*H.*sin(C)).^2));
UP=A.*H.*(-X.^4 - X.^2.*Y.^2.*H.^2 - 2.*X.^2.*F.^2 - X.^2.*H.^2.*cos(C).^2 - 2.*X.*F.*H.^2.*sin(C).*cos(C) - Y.^2.*F.^2.*H.^2 - F.^4 + F.^2.*H.^2.*cos(C).^2 - F.^2.*H.^2).*cos(C)./(H.^2.*(X.*cos(C) + F.*sin(C)).^2 + (X.^2 + Y.^2.*H.^2 + F.^2).^2);

%JJB=JB-UP;

JJB=((X.*H.*cos(C) + F.*H.*sin(C)).*(2.*A.*X.^3 + 2.*A.*X.*F.^2 + 2.*A.*X.*Y.^2.*H.^2 - abs(cos(C)).*imag((X.^4 + F.^4 + 2.*X.^2.*F.^2 - F.^3.*H.*sin(C).*2i - X.^2.*H.^2.*cos(C).^2 - F.^2.*H.^2.*sin(C).^2 - X.^3.*H.*cos(C).*2i - X.*F.^2.*H.*cos(C).*2i - X.^2.*F.*H.*sin(C).*2i - X.*F.*H.^2.*sin(2.*C)).^(1./2)).*abs(A).*abs(H) - A.*X.*H.^2.*cos(C).^2 - (A.*F.*H.^2.*sin(2.*C))./2))./(2.*((F.^2 + Y.^2.*H.^2 + X.^2).^2 + (X.*H.*cos(C) + F.*H.*sin(C)).^2)) - ((F.^2 + Y.^2.*H.^2 + X.^2).*(3.*A.*X.^2.*H.*cos(C) + A.*F.^2.*H.*cos(C) - abs(cos(C)).*real((X.^4 + F.^4 + 2.*X.^2.*F.^2 - F.^3.*H.*sin(C).*2i - X.^2.*H.^2.*cos(C).^2 - F.^2.*H.^2.*sin(C).^2 - X.^3.*H.*cos(C).*2i - X.*F.^2.*H.*cos(C).*2i - X.^2.*F.*H.*sin(C).*2i - X.*F.*H.^2.*sin(2.*C)).^(1./2)).*abs(A).*abs(H) + 2.*A.*X.*F.*H.*sin(C)))./(2.*((F.^2 + Y.^2.*H.^2 + X.^2).^2 + (X.*H.*cos(C) + F.*H.*sin(C)).^2));
surf(X,Y,JJB)
