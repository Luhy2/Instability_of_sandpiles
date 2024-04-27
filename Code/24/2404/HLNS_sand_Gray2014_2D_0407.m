clc
clear

g = 9.8;
sc = deg2rad(29);
F = 1.02;

nu = 1.13e-3;
gamma = 0.517;
beta = 0.136;
CharL = 0.825e-3;
sc1 = deg2rad(20.9);
sc2 = deg2rad(32.76);

hsinf = CharL*gamma*F/beta;
usinf = F*(g*hsinf*cos(sc))^(1/2);

R = usinf*sqrt(hsinf)/nu;
Gamp = gamma*(tan(sc2)-tan(sc1))/(1+gamma)^2;


% N = 211; % -100
N = 1400;
% N = 650;

Nr = N+1;

% --------------------------------
% k2
B = 0;
Minus = -400;
xint = 615-615;
% xend = 1246;
xend = 846-Minus-615;

% N与Nr在上面定义

xx = linspace(xint,xend,Nr);
h = (xend-xint)/(Nr-1);
xreal = xx;

C1 = 34; %1 2 14
C2 = 25;
C3 = 350; % 10

% FA=1+C1*(1+tanh(C2*(-1+(2.*xreal./(xint+xend+C3)))));
FA=ones(1,Nr);
FA=FA';
% FA=1;
% check mesh and FA
% figure()
% plot(xx,'o')
% hold on
% figure()
% plot(xreal,FA,'o')
% hold on

% --------------------------------
% uni dif mat
I = eye(Nr);

fracDl = (-1/(2*h)).*ones(Nr,1);
fracDm = zeros(Nr,1);
fracDr = (1/(2*h)).*ones(Nr,1);

D = spdiags([fracDl fracDm fracDr],-1:1,Nr,Nr);

Dx = full(D);
% forward
Dx(1,1) = -3/(2*h);
Dx(1,2) = 4/(2*h);
Dx(1,3) = -1/(2*h);
% backward
Dx(Nr,Nr) = 3/(2*h);
Dx(Nr,Nr-1) = -4/(2*h);
Dx(Nr,Nr-2) = 1/(2*h);


fracD2l = (1/h^2).*ones(Nr,1);
fracD2m = (-2/h^2).*ones(Nr,1);
fracD2r = (1/h^2).*ones(Nr,1);

D2 = spdiags([fracD2l fracD2m fracD2r],-1:1,Nr,Nr);
D2x = full(D2);
% forward
D2x(1,1) = 2/h^2;
D2x(1,2) = -5/h^2;
D2x(1,3) = 4/h^2;
D2x(1,4) = -1/h^2;
% backward
D2x(Nr,Nr) = 2/h^2;
D2x(Nr,Nr-1) = -5/h^2;
D2x(Nr,Nr-2) = 4/h^2;
D2x(Nr,Nr-3) = -1/h^2;

% --------------------------------
% mat fill
omega = 0.262626;%max growth point
% omega = 0.996315; %中性

LHS11 = Dx-1i*omega*I;
LHS12 = Dx;
LHS21 = -(3/2)*Gamp*I+Dx;
LHS22 = Gamp*I-1i*(F^2)*omega*I+(F^2)*Dx-((F^2)/R)*D2x;

BR = zeros(2*Nr,1);

% BCs

%左边界值H,U,V为设定值
LHS11(1,:) = I(1,:);
% LHS11(1,:) = Dx(1,:);


LHS12(1,:) = 0;

% BR(1,1) = -0.0002271+0.90687i;
BR(1,1) = 0.1072-0.8928i;
% RHS11(1,:) = 0;

LHS22(1,:) = I(1,:);
% LHS22(1,:) = Dx(1,:);


LHS21(1,:) = 0;

% BR(Nr+1,1) = -0.091+0.909i;
BR(Nr+1,1) = 0.1051-0.8753i;
% RHS22(1,:) = 0;


% 右边界无条件
% LHS21(Nr,:) = 0;
% LHS22(Nr,:) = I(Nr,:);
% BR(2*Nr-1,1) = 0;
% BR(2*Nr,1) = 0;

LHS = [LHS11,LHS12;LHS21,LHS22];

x1 = linsolve(LHS,BR);


Hreal = real(x1(1:Nr,1));
Himag = imag(x1(1:Nr,1));


% figure()
subplot(3,1,1)
plot(xreal,real(x1(1:Nr,1)),'-o')
hold on
plot(xreal,imag(x1(1:Nr,1)),'-')
hold on
plot(xreal,abs(x1(1:Nr,1)),'--','LineWidth',3)
subplot(3,1,2)
plot(xreal,real(x1(Nr+1:2*Nr,1)),'-o')
hold on
plot(xreal,imag(x1(Nr+1:2*Nr,1)),'-')
hold on
plot(xreal,abs(x1(1:Nr,1)),'--','LineWidth',3)
subplot(3,1,3)
plot(xreal,FA,'-o')
hold on
