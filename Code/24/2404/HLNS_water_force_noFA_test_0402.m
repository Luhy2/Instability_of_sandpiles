clc
clear



g=9.8;
sc=deg2rad(4.6);
rint=2.5e-3;

%q0=3e-5/(2*pi);
q0 = 3.534e-3/(2*pi);

% aint=80*rint;
% aend=2000*rint;
aint = 80*rint;
aend = 2000*rint;

c0=q0;

nu=4.89e-6;

C=c0;
c1=27/35;
c2=3;

EXT=c2*nu*C;
hint=(C*c2*nu/(g*sin(sc)*aend))^(1/3);


N = 205;
% N = 650;

Nr = N+1;

a = linspace(aint,aend,Nr);
hs = (C*c2*nu./(g*sin(sc).*a)).^(1/3);
us = C./(a.*hs);
% test base flow
% plot(a,hs,'o')
% hold on
% plot(a,us,'o')

usinf = us(1,end);
hsinf = hs(1,end);

R = usinf*hsinf/nu;
F = usinf/(g*hsinf*cos(sc))^(1/2);
% --------------------------------
% k2
B = 0.9;

xint = 615;
% xend = 1246;
xend = 846;

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

Spterm = 1i*B./(xx.*cos(sc));
% Spterm = 1*B./(xx.*cos(sc));
SPT = Spterm.*I;

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
% wf = 0.4716;
wf = 0.911414+0.0308013i;

LHS11 = Dx-1i*wf*I;
LHS12 = Dx;
LHS13 = SPT;
LHS21 = (1/F^2)*Dx-(6/R).*I;
LHS22 = Dx-(1/R)*D2x.*FA+(3/R).*I-1i*wf*I;
LHS23 = zeros(Nr,Nr);
LHS31 = (1/F^2)*SPT;
LHS32 = zeros(Nr,Nr);
LHS33 = Dx-(1/R)*D2x.*FA+(3/R).*I-1i*wf*I;

% RHS11 = I;
% RHS12 = zeros(Nr,Nr);
% RHS13 = zeros(Nr,Nr);
% RHS21 = zeros(Nr,Nr);
% RHS22 = I;
% RHS23 = zeros(Nr,Nr);
% RHS31 = zeros(Nr,Nr);
% RHS32 = zeros(Nr,Nr);
% RHS33 = I;

BR = zeros(3*Nr,1);

% BCs

%左边界值H,U,V为设定值
LHS11(1,:) = I(1,:);
% LHS11(1,:) = Dx(1,:);
% LHS11(1,:) = Dx(1,:)-1i*0.407*I(1,:);

LHS12(1,:) = 0;
LHS13(1,:) = 0;
% BR(1,1) = 0.6109-0.0594i; 
BR(1,1) = 0.08;
% RHS11(1,:) = 0;

LHS22(1,:) = I(1,:);
% LHS22(1,:) = Dx(1,:);
% LHS22(1,:) = Dx(1,:)-1i*0.407*I(1,:);

LHS21(1,:) = 0;
LHS23(1,:) = 0;
% BR(2,1) = 0.7895;
BR(Nr+1,1) = 0.1;
% RHS22(1,:) = 0;

LHS33(1,:) = I(1,:);
% LHS33(1,:) = Dx(1,:);
% LHS33(1,:) = Dx(1,:)-1i*0.407*I(1,:);

LHS31(1,:) = 0;
LHS32(1,:) = 0;
BR(3,1) = 0;
% RHS33(1,:) = 0;

% 右边界值U,V为0
% LHS11(Nr,:) = I(Nr,:);
% LHS12(Nr,:) = 0;
% RHS11(Nr,:) = 0;

% LHS22(Nr,:) = I(Nr,:);
% LHS22(Nr,:) = Dx(Nr,:);
% LHS21(Nr,:) = 0;
% LHS23(Nr,:) = 0;
% RHS22(Nr,:) = 0;
BR(3*Nr-1,1) = 0;

% LHS33(Nr,:) = I(Nr,:);
% LHS33(Nr,:) = Dx(Nr,:);
% LHS31(Nr,:) = 0;
% LHS32(Nr,:) = 0;
% RHS33(Nr,:) = 0;
BR(3*Nr,1) = 0;

LHS = [LHS11,LHS12,LHS13;LHS21,LHS22,LHS23;LHS31,LHS32,LHS33];
% RHS = 1i.*[RHS11,RHS12,RHS13;RHS21,RHS22,RHS23;RHS31,RHS32,RHS33];

x1 = linsolve(LHS,BR);
% figure()
subplot(4,1,1)
plot(xreal,real(x1(1:Nr,1)),'-o')
hold on
subplot(4,1,2)
plot(xreal,real(x1(Nr+1:2*Nr,1)),'-o')
hold on
subplot(4,1,3)
plot(xreal,real(x1(2*Nr+1:end,1)),'-o')
hold on
subplot(4,1,4)
plot(xreal,FA,'o')
hold on
