clc
clearvars -except what_you_want_to_remain
% define parameters in the disturbance case (V1.0)

% I=i;  %legal symbol
% A=k_1;
% B=k_2;
% C=\omega;
% D=u_s;
% E=h_s;
% F=\cos{\alpha};
% G=g;
% H=\frac{1}{h_2}; %h_2=a*\cos{\alpha}
% P=p; %p=(\frac{\partial Y}{\partial u})_s
% Q=q; %q=(\frac{\partial Y}{\partial_h}+p*\frac{\partial_u}{\partial_h})_s
% Y=Y_s; %Constitutive relation and velocity profile determined term

%A,B,C,D,E,F,G,H = syms('A,B,C,D,E,F,G,H','real');
%syms A B C D E F G H real
syms  A B D E F G H P Q Y real

%X = sym('X');
%A=sym('A');
%B=sym('B');
C=sym('C');
    %test
    %UPI=A+B*C*sin(D);
    %UPI
I=1i;

CW=54/35;
P=0;

M = [I*A*D-I*C,I*A*E,I*B*E*H;
    I*A*G*F-Q*D,I*A*D-I*C-P*D-Y,0;
    I*B*H*G*F,0,I*A*D-I*C-Y
];
%latex(M)
KK=det(M);
latex(KK);
HS=solve(KK,C)
latex(HS);

%RESULT1=simplify(HS(1,1));
%RESULT2=simplify(HS(2,1));

%RE1=real(RESULT1)
%IM1=imag(RESULT1)
%RE2=real(RESULT2)
%IM2=imag(RESULT2)
%%
clc
clear
nu=50e-6;%kinematic viscosity \nu=\mu/\rho
g=9.8;
sc=deg2rad(30);
rint=2.5e-3;
q0=3e-5/(2*pi);
aint=100*rint;

hs=((3*nu*q0/(g*sin(sc)))*(1/aint))^(1/3);
us=(1/3)*(g*sin(sc)/nu)*(hs)^2;
Ys=3*nu/(hs)^2;
h2=aint*cos(sc);

D=us;
E=hs;
F=cos(sc);
G=g;
H=1/h2;
P=0;
Q=-6*nu/(hs)^3;
Y=Ys;


%syms  A B real
A=1;
B=1;


C=sym('C');

I=1i;

M = [I*A*D-I*C,I*A*E,I*B*E*H;
    I*A*G*F-Q*D,I*A*D-I*C-P*D-Y,0;
    I*B*H*G*F,0,I*A*D-I*C-Y
];
%latex(M)

KK=det(M);
latex(KK);
%%KK=eig(M)


HS=solve(KK,C);
latex(HS);
IMHS=imag(HS);
latex(IMHS)

   


