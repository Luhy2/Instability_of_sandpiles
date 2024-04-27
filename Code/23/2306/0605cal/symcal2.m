clc
clearvars -except QWEHG
% define parameters (V4.0)

% L=\lambda 关键字使用其他字母做替换
% I=i 虚数单位单独替换
% A=k_a
% Q=u_a^*
% W=\partial_a{u_a^*}
% E=1/a
% H=h^*
% T=k_{\theta}
% G=g\cos{\alpha}

%A,B,C,D,E,F,G,H = syms('A,B,C,D,E,F,G,H','real');
%syms A B C D E F G H real
syms  A Q W E H T G real

%X = sym('X');
L=sym('L');
    %test
    %UPI=A+B*C*sin(D);
    %UPI
    
I=1i;

%M =[X+I*D*A+H*I*E*G,-2*H*cos(C)*G,0,B*I*D ; cos(C)*H*G,X+I*D*A+H*I*E*G+cos(C)*H*A,sin(C)*H*G,B*H*I*E ; 0,-2*H*cos(C)*G,X+I*D*A+H*I*E*G,B*I*F ; I*D+cos(C)*H,I*H*E,I*F+sin(C)*H,0];
M = [L+I*A*Q+W+Q*E,I*A*H,I*T*H*E;
    I*A*G,L+I*A*Q+W,0;
    I*T*E*G,0,L+I*A*Q+Q*E];

KK=det(M);

%latex(KK)
%Q=11.00358;
%W=0.409;
%E=0.315671;
%H=0.001434;
%G=8.487;


HS=solve(KK,L)

%RESULT1=simplify(HS(1,1));
%RESULT2=simplify(HS(2,1));

%RE1=real(RESULT1)
%IM1=imag(RESULT1)
%RE2=real(RESULT2)
%IM2=imag(RESULT2)
