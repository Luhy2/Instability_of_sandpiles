%%
% first test dH+U+H=0 ; d2U+U=0
clc
clear

xint = 0;
%xend = 2*pi;% real range 0 - 2pi(6.28)
xend = 200*pi;

N = 1250;
Nr = N+1;

I = eye(Nr);

xx = linspace(xint,xend,Nr);
h = (xend-xint)/(Nr-1);
xreal = xx;

% FA settings
C1 = 1; %34 basic
C2 = 7;
C3 = 50; % 10

FA2=1+C1*(1+tanh(C2*(-1+(2.*xx./(xint+xend+C3)))));
% FA1 = 1+C1*(1+tanh(C2*(-1+(2.*flip(xx1)./(x1int+xint+C3)))));
% FA = [FA1(1:end-1),FA2];
FA = FA2;
plot(xreal,FA2,'o')
hold on
%%
% uni dif mat

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

% xreal = (xint+xend)/2 + xx*(xend-xint)/2;
% xreal = flip(xreal);

% check mesh
% figure()
% plot(xx,'o')
% hold on
% plot(xreal,'^')

% mat fill

LHS11 = Dx;
LHS12 = I;
LHS21 = zeros(Nr,Nr);
LHS22 = FA.*D2x;

RHS11 = -I;
% RHS11 = zeros(Nr,Nr);
RHS12 = zeros(Nr,Nr);
RHS21 = zeros(Nr,Nr);
RHS22 = -I;

% BCs

%左边界值为0
LHS11(1,:) = I(1,:);
LHS12(1,:) = 0;
RHS11(1,:) = 0;

LHS22(1,:) = I(1,:);
LHS21(1,:) = 0;
RHS22(1,:) = 0;

% 右边界值为0
% LHS11(Nr,:) = I(Nr,:);
% LHS12(Nr,:) = 0;
% RHS11(Nr,:) = 0;

LHS22(Nr,:) = I(Nr,:);
LHS21(Nr,:) = 0;
RHS22(Nr,:) = 0;

LHS = [LHS11,LHS12;LHS21,LHS22];
RHS = [RHS11,RHS12;RHS21,RHS22];
%%
% solve
[V,D] = eig(LHS,RHS);
% [V,D] = eigs(LHS,RHS,20,'smallestabs');

Omegalist=diag(D);
Omeganew=Omegalist;
Vnew=V;


realOme=real(Omeganew);
imagOme=imag(Omeganew);

figure()
plot(realOme,imagOme,'o')
xlim([0,10e-4])
hold on
for i = 1:1:length(Omeganew)
    text(realOme(i), imagOme(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end
%%

TRY = 2496;

C=realOme(TRY,1);
%Aal = 0.0097;%2498
Aal = 0.0037;
hreal = Aal.*(-exp(-C.*xx) - sqrt(C).*sin(sqrt(C).*xx) + cos(sqrt(C).*xx))./(sqrt(C).*(C + 1));
ureal = Aal.*sin(sqrt(C).*xx);


figure()
subplot(2,1,1)
plot(xreal,real(Vnew(1:Nr,TRY)),'-o')
hold on
% ylim([0,2*10e-13])
plot(xreal,imag(Vnew(1:Nr,TRY)),'-')
hold on
plot(xx,-hreal,'m-','LineWidth',2)
hold on

subplot(2,1,2)
plot(xreal,real(Vnew(Nr+1:end,TRY)),'-o')
hold on
plot(xreal,imag(Vnew(Nr+1:end,TRY)),'-')
hold on

plot(xx,-ureal,'m-','LineWidth',2)
hold on






