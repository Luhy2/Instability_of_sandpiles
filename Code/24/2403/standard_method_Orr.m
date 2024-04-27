clc;
clear all;
nosmod = 300;
Re = 5773;
alpha = 1.02;
beta = 0;
iflow = 1;

iUnit=1i;

[D,x]=cheb0(nosmod);
D1=D;
D2=D1^2;
D4=D2^2;

nMatrix=nosmod+1;

k=sqrt(alpha^2+beta^2);
I = eye(nMatrix);

% -----------------------------
% velocity of Poiseuille
% -----------------------------
U = 1-x.^2;
UDiag=diag(U);
Udiff2=-2*I;
Udiff1=-2*x;
Udiff1Diag=diag(Udiff1);
Udiff2Diag=(Udiff2);

% set up OS and SQ matrix
nabla2=D2-alpha^2*I;
%nabla4=D4-2*D2*I*k^2+alpha^2*I;
nabla4=D4-2*k^2*D2+k^4*I;
% -----------------------------
% B
% -----------------------------
B11 = [nabla2];
B12= zeros(nMatrix,nMatrix);
B21= zeros(nMatrix,nMatrix);
B22 = -I;

B11(1,:)=0;
B11(2,:)=0;
B11(nMatrix-1,:)=0;
B11(nMatrix,:)=0;
B22(1,:)=0;
B22(nMatrix,:)=0;
B=[B11,B12;B21,B22];

% -----------------------------
% A
% -----------------------------
LOS = 1/Re*nabla4   + iUnit*alpha*Udiff2Diag*I  - iUnit*alpha*UDiag * nabla2;
LSQ = iUnit*alpha*UDiag-1/Re*nabla2;
A11 = LOS;
A12 = zeros(nMatrix,nMatrix);
A21 = iUnit*beta*Udiff1Diag;
A22 = LSQ;

para=1;
A11(1,:)=para*I(1,:);
A11(2,:)=para*D1(1,:);
A11(nMatrix-1,:)=para*D1(nMatrix,:);
A11(nMatrix,:)=para*I(nMatrix,:);
A21(1,:)=0;
A21(nMatrix,:)=0;
A22(1,:)=para*I(1,:);
A22(nMatrix,:)=para*I(nMatrix,:);
A=[A11,A12;A21,A22];

% -----------------------------
% Solving eigenvalue problems
% -----------------------------
% lam = eig(A11,B11);
[V,lam] = eig(A,B);

lam =iUnit* diag(lam);

RL = real(lam(1:2*(nMatrix-1),1));
IL = imag(lam(1:2*(nMatrix-1),1));

plot(RL,IL,'o')
hold on
% for i = 1:1:length(lam)-2
%     text(RL(i), IL(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
% end

xlim([0,1])
ylim([-1,0.1])
%%
figure()
TRY = 365;
plot(real(V(:,TRY)),'-o')
hold on
plot(imag(V(:,TRY)),'-o')
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

