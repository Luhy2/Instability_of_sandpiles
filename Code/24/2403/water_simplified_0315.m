clc
clear

% init para
k2 = 2.6;

g = 9.8;
sc = deg2rad(20);
rint = 2.5e-3;
Qstar = 3e-5; % init Q flow rate
qstar = Qstar/(2*pi);
uint = (2*qstar)/(rint^2);
q0 = qstar/(uint*(rint^2));

nu = 50e-6; % kinematic viscosity \nu=\mu/\rho

R0 = uint*rint/nu; % 2*R0 = Rej
F0 = sqrt((uint^2)/(g*rint*cos(sc)));

c1 = 2*27/35;
c2 = 3;

    % create standard base flow H(A) and U(A) distribution
    Aint = 280;
    Aend = 20;
    Nspace = 200;
    AA = linspace(Aint,Aend,Nspace);
    
    HA = (((F0^2)/R0)*c2*q0./(AA.*tan(sc))).^(1/3);
    UUA = q0./(HA.*AA);
    
% create mesh
[D,x] = cheb0(Nspace);
D2 = D^2;
x = x';
    % remap
    Xspace =((Aend+Aint)/2)+((Aint-Aend)/2).*x;
    Hcheb = (((F0^2)/R0)*c2*q0./(Xspace.*tan(sc))).^(1/3);
    Ucheb = q0./(Hcheb.*Xspace);

% plot(flip(AA),'o')
% hold on
% plot(Xspace,'o')
% hold on

%     subplot(2,1,1)
%     plot(AA,HA,'^')
%     hold on
%     plot(Xspace,Hcheb,'o')
%     hold on
%     subplot(2,1,2)
%     plot(AA,UUA,'^')
%     hold on
%     plot(Xspace,Ucheb,'o')
%     hold on

% figure()
% plot(x,Hcheb)

% cancel first and last row and col
xcan = x(1,2:Nspace);
Xspacecan = Xspace(1,2:Nspace);
D = D(2:Nspace,2:Nspace);
D2 = D2(2:Nspace,2:Nspace);
    Hcheb = (((F0^2)/R0)*c2*q0./(Xspacecan.*tan(sc))).^(1/3);
    Ucheb = q0./(Hcheb.*Xspacecan);


LHS11 = diag(Ucheb)*D;
[Nrow,Ncol] = size(LHS11);
% OneM = ones(Nrow,Ncol);
EYEM = eye(Nrow);

LHS12 = diag(Hcheb)*D;
LHS13 = 1i*k2*diag(Hcheb);
% LHS21 = (1/F0^2)*D-(c2/R0).*(1./(diag(Hcheb)).^2).*(2.*diag(Ucheb)./diag(Hcheb));
LHS21 = (1/F0^2)*D-(2*c2/R0)*diag(Ucheb./(Hcheb.^3));
LHS22 = diag(Ucheb)*D+(1/R0)*c2*diag(1./(Hcheb.^2))+(1/R0)*D2;
LHS23 = 0*EYEM;
LHS31 = (1/F0^2)*1i*k2*EYEM;
LHS32 = 0*EYEM;
LHS33 = diag(Ucheb)*D+(1/R0)*c2*diag(1./(Hcheb.^2))+(1/R0)*D2;

MATSIZE = 3*(Nspace-1);
LLHS = zeros(MATSIZE);

for i = 1:Nrow
    for j = 1:Ncol
        LLHS(1+3*(i-1),1+3*(j-1)) = LHS11(i,j);
        LLHS(1+3*(i-1),2+3*(j-1)) = LHS12(i,j);
        LLHS(1+3*(i-1),3+3*(j-1)) = LHS13(i,j);
        LLHS(2+3*(i-1),1+3*(j-1)) = LHS21(i,j);
        LLHS(2+3*(i-1),2+3*(j-1)) = LHS22(i,j);
        LLHS(2+3*(i-1),3+3*(j-1)) = LHS23(i,j);
        LLHS(3+3*(i-1),1+3*(j-1)) = LHS31(i,j);
        LLHS(3+3*(i-1),2+3*(j-1)) = LHS32(i,j);
        LLHS(3+3*(i-1),3+3*(j-1)) = LHS33(i,j);
    end
end

RRHS = eye(3*Nrow);

    
% [V,D] = eig(LLHS,RRHS);
[V,Dz] = eigs(LLHS,RRHS,500,'smallestabs');
Omegalist = diag(Dz)./1i;
RL = real(Omegalist);
IM = imag(Omegalist);

%%
plot(RL,IM,'o')
hold on
% for i = 1:1:length(Omegalist)
%     text(RL(i), IM(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
% end
%%
figure()
clc
TRY=1;
for i=1:(MATSIZE/3)
    AlphaH(i,1)=V((i-1)*3+1,TRY);
    AlphaU(i,1)=V((i-1)*3+2,TRY);
    AlphaV(i,1)=V((i-1)*3+3,TRY);
end

subplot(3,1,1)
plot(x(2:end-1),real(AlphaH),'-o')
hold on
plot(x(2:end-1),imag(AlphaH))
hold on

subplot(3,1,2)
plot(x(2:end-1),real(AlphaU),'-o')
hold on
plot(x(2:end-1),imag(AlphaU))
hold on

subplot(3,1,3)
plot(x(2:end-1),real(AlphaV),'-o')
hold on
plot(x(2:end-1),imag(AlphaV))
hold on

%%

% CHEB compute D = differentiation matrix, x = Chebyshev grid
function [ D,x ] = cheb0( N )

if N==0, 
    D = 0;
    x=1;
    return,
end
x = -cos(pi * (0:N)/N)';
c = [2; ones(N-1,1);2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X - X';
D = (c*(1./c)')./(dX + (eye(N+1)));          % off-diagonal entries
D = D - diag(sum(D'));                       % diagonal entries
end

