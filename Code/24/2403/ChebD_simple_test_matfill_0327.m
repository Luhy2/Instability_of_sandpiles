clc
clear

% y''=-ky reduction mat fill test chebD
xint = 0;
xm1 = 20;

N = 50;
[D,x] = cheb(N);

Dx = D/((xm1-xint)/2);
% Dx = D;

Dx2 = Dx^2;

Nr = N+1;
I = eye(Nr);

% uniform mesh
xx = linspace(xint,xm1,Nr);

% mesh remap
xreal = (xint+xm1)/2 + x*(xm1-xint)/2;
xreal = flip(xreal);

C1 = 7; %1 2 14
C2 = 15;
C3 = 15; % 10

FA=1+C1*(1+tanh(C2*(-1+(2.*xreal./(xint+xm1+C3)))));
% FA = 1;

% examine mesh and FA

%     plot(xreal,'o')
%     hold on
%     plot(xx)
%     hold on
% plot(xreal,FA,'-o')
% hold on

% LHS mat fill

LHS11 = zeros(Nr,Nr);
LHS12 = FA.*Dx;
LHS21 = Dx;
LHS22 = -I;

RHS11 = -I;
RHS12 = zeros(Nr,Nr);
RHS21 = zeros(Nr,Nr);
RHS22 = zeros(Nr,Nr);

%边界条件设置
%左边界值为0
LHS11(1,:) = I(1,:);
LHS12(1,:) = 0;
RHS11(1,:) = 0;

% 右边界一阶导为0
% LHS21(Nr,:) = 0;
% LHS22(Nr,:) = I(Nr,:);


% 右边界值为0
% 填在上矩阵的方式
LHS11(Nr,:) = I(Nr,:);
LHS12(Nr,:) = 0;
RHS11(Nr,:) = 0;

%填在下矩阵的方式 不好
% LHS21(Nr,:) = 0;
% LHS22(Nr,:) = D(Nr,:);

    
LHS = [LHS11,LHS12;LHS21,LHS22];
RHS = [RHS11,RHS12;RHS21,RHS22];

[V,D] = eigs(LHS,RHS,20,'smallestabs');
% [V,D] = eig(LHS,RHS);

Omegalist=diag(D);
Omeganew=Omegalist;
Vnew=V;

figure()
realOme=real(Omeganew);
imagOme=imag(Omeganew);
% plot(realOme,imagOme,'o')
% plot(realOme,'o')
% hold on
% for i = 1:1:length(Omeganew)
%     text(realOme(i), imagOme(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
% end

for i = 1:1:length(Omeganew)
    text(i,realOme(i),num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end


TRY = 4;

% figure()
plot(xreal,Vnew(1:Nr,TRY),'-o')
hold on

% figure()
% plot(xreal,Vnew(Nr+1:end,TRY),'-o')
% hold on




