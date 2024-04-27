clc
clear


LHS = [2,1,-5,1;
    1,-5,0,7;
    0,2,1,-1;
    1,6,-1,-4];

% RHS = [13,-9,6,0];
RHS = [0,0,0,0];
RHS = RHS';

x1 = LHS\RHS;

x2 = linsolve(LHS,RHS);
x1
x2
res = x1-x2;
res

%%
clc
IMAG = [1i,0;
    0,1i]

IMAGinv1 = IMAG^(-1)
IMAGinv2 = inv(IMAG)

%%
clc
sc = deg2rad(4.6);
sc

sin(sc)

%%
clc
nu = 4.89e-6;
nuu = 4.89*10^(-6);

res = nuu-nu
%%
clc
MAT = [2,-1;1,4];
% MAT = [-2,1;-1,-4];
[V,D]=eig(MAT)

%%
clc
clear
RHS = [0,0]';
LHS = [-1,-1;1,1];

x1 = linsolve(LHS,RHS)

%%
OP = [2,1;1,2];

2*OP

2.*OP-1*eye(2)
%%
clc
clear

A = [1,2;3,4];
B = 2*A;
BB = 2.*A;

Two = 2.*eye(2);
B3 = Two*A;
B4 = Two.*A;





















