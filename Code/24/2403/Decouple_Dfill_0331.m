clc
clear

xint = 0;
xend = 2*pi;

N = 250;
Nr = N+1;

[D,x] = cheb(N);

Dx = D/((xend-xint)/2);
I = eye(Nr);

Dx2 = Dx^2;

xx = linspace(xint,xend,Nr);

xreal = (xint+xend)/2 + x*(xend-xint)/2;
xreal = flip(xreal);

% check mesh
% plot(xx,'o')
% hold on
% plot(xreal,'^')

% mat fill

LHS11 = Dx2;
LHS12 = zeros(Nr,Nr);
LHS21 = I;
LHS22 = Dx2;

RHS11 = -I;
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
LHS11(Nr,:) = I(Nr,:);
LHS12(Nr,:) = 0;
RHS11(Nr,:) = 0;

LHS22(Nr,:) = I(Nr,:);
LHS21(Nr,:) = 0;
RHS22(Nr,:) = 0;

LHS = [LHS11,LHS12;LHS21,LHS22];
RHS = [RHS11,RHS12;RHS21,RHS22];

% solve
[V,D] = eig(LHS,RHS);
% [V,D] = eigs(LHS,RHS,20,'smallestabs');

Omegalist=diag(D);
Omeganew=Omegalist;
Vnew=V;


realOme=real(Omeganew);
imagOme=imag(Omeganew);
%%
figure()
plot(realOme,imagOme,'o')
xlim([-10,30])
hold on
for i = 1:1:length(Omeganew)
    text(realOme(i), imagOme(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end
%%

TRY = 449;

figure()
subplot(2,1,1)
plot(xreal,real(Vnew(1:Nr,TRY)),'-o')
hold on
plot(xreal,imag(Vnew(1:Nr,TRY)),'-')
hold on
subplot(2,1,2)
plot(xreal,real(Vnew(Nr+1:end,TRY)),'-o')
hold on
plot(xreal,imag(Vnew(Nr+1:end,TRY)),'-')
hold on

