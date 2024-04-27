clc
clear
% define parameters (2)

% A=u_a^*=u_1
% B=\frac{1}{\rho\phi}
% C=\alpha
% D=k_a
% E=k_{\theta}
% F=k_b

% G=u_{\theta}^*=u_2

% H=\frac{1}{h_2}
% X=\lambda 关键字使用其他字母做替换
%A,B,C,D,E,F,G,H = syms('A,B,C,D,E,F,G,H','real');
syms A B C D E F G H real

X = sym('X');

    %test
    %UPI=A+B*C*sin(D);
    %UPI
    
I=1i;

M =[X+I*D*A+H*I*E*G,-2*H*cos(C)*G,0,B*I*D ; cos(C)*H*G,X+I*D*A+H*I*E*G+cos(C)*H*A,sin(C)*H*G,B*H*I*E ; 0,-2*H*cos(C)*G,X+I*D*A+H*I*E*G,B*I*F ; I*D+cos(C)*H,I*H*E,I*F+sin(C)*H,0];


KK=det(M);

HS=solve(KK,X);


RESULT1=simplify(HS(1,1));
RESULT2=simplify(HS(2,1));

%RE1=real(RESULT1)
%IM1=imag(RESULT1)
%RE2=real(RESULT2)
%IM2=imag(RESULT2)



%% benchmark
%GHS=subs(HS,G,0);

%SGHS=simplify(GHS);
%latex(SGHS)

%RIO=SGHS(1,1);
%JB=real(RIO);
%latex(JB)



