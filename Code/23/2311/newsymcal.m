clc
clearvars -except what_you_want_to_remain
% define parameters in the disturbance case (V2.0)

% I=i;  %legal symbol
% A=k_1;
% B=k_2;
% C=\omega;
% T=5/4; %profile coeff

% F=(1/Fr^2);
% M=3/2;
% H=\frac{1}{X_2}; %h_2=a*\cos{\alpha}%X_2=A*\cos{\alpha}
% Y=\Gamma; %bed stress coeff
% R=(1/Re*)*(k1^2+(1/x2)^2 *k2^2) %effective Reynold number



syms  A B F H T M R Y real

C=sym('C');
    
I=1i;

M = [I*A-I*C,I*A,I*B*H;
    F*(I*A-M*Y),T*I*A-I*C+R+F*Y,I*B*H;
    F*(I*B*H-M*Y),0,I*A-I*C+R+F*Y
];
%latex(M)
KK=det(M);
%latex(KK)
HS=solve(KK,C);
latex(HS)

%RESULT1=simplify(HS(1,1));
%RESULT2=simplify(HS(2,1));

%RE1=real(RESULT1)
%IM1=imag(RESULT1)
%RE2=real(RESULT2)
%IM2=imag(RESULT2)
