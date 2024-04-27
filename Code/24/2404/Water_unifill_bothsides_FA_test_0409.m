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


% N = 449;
N = 649;
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
B = 0;

% x1int = 384;
x1int = 5;
xint = 615;
% xend = 1246;
% xend = 846;
xend = 615+610;

% N与Nr在上面定义

xx = linspace(xint,xend,Nr/2);
h = (xend-xint)/(Nr/2-1);

xx1 = linspace(x1int-h,xint,Nr/2+1);
h1 = (xint-x1int+h)/(Nr/2+1-1);

xreal = [xx1(1:end-1),xx];
%xreal = xx;

C1 = 10; %34 basic
C2 = 45;
C3 = 450; % 10

FA2=1+C1*(1+tanh(C2*(-1+(2.*xx./(xint+xend+C3)))));
FA1 = 1+C1*(1+tanh(C2*(-1+(2.*flip(xx1)./(x1int+xint+C3)))));
FA1 = ones(1,Nr/2+1);
FA = [FA1(1:end-1),FA2];
% FA=ones(1,Nr);
FA=FA';
% FA=1;
% check mesh and FA
% figure()
% plot(xx,'o')
% hold on
figure()
plot(xreal,FA,'o')
hold on

% --------------------------------
% uni dif mat
I = eye(Nr);

Spterm = 1i*B./(xreal.*cos(sc));
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

LHS11 = Dx;
LHS12 = Dx;
LHS13 = SPT;
LHS21 = (1/F^2)*Dx-(6/R).*I.*FA;
LHS22 = Dx-(1/R)*D2x.*FA+(3/R).*I.*FA;
LHS23 = zeros(Nr,Nr);
LHS31 = (1/F^2)*SPT;
LHS32 = zeros(Nr,Nr);
LHS33 = Dx-(1/R)*D2x.*FA+(3/R).*I.*FA;

RHS11 = I;
RHS12 = zeros(Nr,Nr);
RHS13 = zeros(Nr,Nr);
RHS21 = zeros(Nr,Nr);
RHS22 = I;
RHS23 = zeros(Nr,Nr);
RHS31 = zeros(Nr,Nr);
RHS32 = zeros(Nr,Nr);
RHS33 = I;

% BCs
% wf = 0.407;
%左边界值H,U,V为0
LHS11(1,:) = I(1,:);
LHS11(1,Nr) = -1;
% LHS11(1,:) = Dx(1,:);
% LHS11(1,:) = Dx(1,:)-1i*wf*I(1,:);

LHS12(1,:) = 0;
LHS13(1,:) = 0;
RHS11(1,:) = 0;

LHS22(1,:) = I(1,:);
LHS22(1,Nr) = -1;
% LHS22(1,:) = Dx(1,:);
% LHS22(1,:) = Dx(1,:)-1i*wf*I(1,:);

LHS21(1,:) = 0;
LHS23(1,:) = 0;
RHS22(1,:) = 0;

LHS33(1,:) = I(1,:);
LHS33(1,Nr) = -1;
% LHS33(1,:) = Dx(1,:);
% LHS33(1,:) = Dx(1,:)-1i*wf*I(1,:); %0.207

LHS31(1,:) = 0;
LHS32(1,:) = 0;
RHS33(1,:) = 0;

% 右边界值U,V为0
% LHS11(Nr,:) = I(Nr,:);
% LHS11(Nr,1) = -1;
% LHS12(Nr,:) = 0;
% LHS13(Nr,:) = 0;
% RHS11(Nr,:) = 0;

LHS22(Nr,:) = I(Nr,:);
LHS22(Nr,1) = -1;
% LHS22(Nr,:) = Dx(Nr,:);
LHS21(Nr,:) = 0;
LHS23(Nr,:) = 0;
RHS22(Nr,:) = 0;

LHS33(Nr,:) = I(Nr,:);
LHS33(Nr,1) = -1;
% LHS33(Nr,:) = Dx(Nr,:);
LHS31(Nr,:) = 0;
LHS32(Nr,:) = 0;
RHS33(Nr,:) = 0;

LHS = [LHS11,LHS12,LHS13;LHS21,LHS22,LHS23;LHS31,LHS32,LHS33];
RHS = 1i.*[RHS11,RHS12,RHS13;RHS21,RHS22,RHS23;RHS31,RHS32,RHS33];
%%
% --------------------------------
% cal eig and vectors

[V,Dnew] = eig(LHS,RHS);
% [V,D] = eigs(LHS,RHS,60,'smallestabs');

Omegalist=diag(Dnew);
Omeganew=Omegalist;
Vnew=V;

realOme=real(Omeganew);
imagOme=imag(Omeganew);

figure()
plot(realOme,imagOme,'o')
% ylim([-10e-3,0])
hold on
for i = 1:1:length(Omeganew)
    text(realOme(i), imagOme(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end
%%
TRY = 244;
figure()
plot(xreal,real(Vnew(Nr+1:2*Nr,TRY)),'-o')
hold on

%%
TRY = 406;
figure()
subplot(3,1,1)
plot(xreal,real(Vnew(1:Nr,TRY)),'-o')
hold on
plot(xreal,imag(Vnew(1:Nr,TRY)),'-')
hold on
% plot(xreal,abs(Vnew(1:Nr,TRY)),'-^')
% hold on

subplot(3,1,2)
plot(xreal,real(Vnew(Nr+1:2*Nr,TRY)),'-o')
hold on
plot(xreal,imag(Vnew(Nr+1:2*Nr,TRY)),'-')
hold on

subplot(3,1,3)
% plot(xreal,real(Vnew(2*Nr+1:end,TRY)),'-o')
% hold on
% plot(xreal,imag(Vnew(2*Nr+1:end,TRY)),'-')
% hold on
plot(xreal,FA,'o')
hold on

% Output(:,2) = real(Vnew(1:Nr,TRY))

