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
    % Last time R denoted as R0
% R=1/Re*
    % R0=R*(A^2+(H^2)*(B^2))

%syms  A B F H T M R Y real
%C=sym('C');

syms B C F H T M R Y real
A=sym('A');

I=1i;

%M = [I*A-I*C,I*A,I*B*H;
%    F*(I*A-M*Y),T*I*A-I*C+R+F*Y,I*B*H;
%    F*(I*B*H-M*Y),0,I*A-I*C+R+F*Y
%];
M = [I*A-I*C,I*A,I*B*H;
    F*(I*A-M*Y),T*I*A-I*C+R*(A^2+(H^2)*(B^2))+F*Y,I*B*H;
    F*(I*B*H-M*Y),0,I*A-I*C+R*(A^2+(H^2)*(B^2))+F*Y
];
%latex(M)
KK=det(M);
%latex(KK)
%HS=solve(KK,C);
HS=solve(KK,A);

latex(HS)

%re=simplify(HS);
%latex(re)
%RESULT1=simplify(HS(1,1));
%RESULT2=simplify(HS(2,1));

%RE1=real(RESULT1)
%IM1=imag(RESULT1)
%RE2=real(RESULT2)
%IM2=imag(RESULT2)

%%
syms R T F C M B H z

% Define the polynomial expression
polynomial = R^2*z^5 + R*T*z^4*1i - F*R*z^4*1i - C*R^2*z^4 + R*z^4*1i + F*M*R*Y*z^3 + 2*B^2*H^2*R^2*z^3 + 2*F*R*Y*z^3 - C*R*T*z^3*1i - C*R*z^3*3i - T*z^3 + F*z^3 + B*F*H*M*R*Y*z^2 - 2*B^2*C*H^2*R^2*z^2 - 2*C*F*R*Y*z^2 + B^2*H^2*R*T*z^2*1i - B^2*F*H^2*R*z^2*2i + B^2*H^2*R*z^2*1i + F*T*Y*z^2*1i + F*M*Y*z^2*1i + F*Y*z^2*1i + 2*C*T*z^2 - C*F*z^2 - F^2*Y*z^2*1i + C^2*R*z^2*2i + C*z^2 + B*F*H*M*T*Y*z*1i + B^2*F*H^2*M*R*Y*z - B*F*H*M*Y*z*1i + 2*B^2*F*H^2*R*Y*z - B^2*C*H^2*R*T*z*1i - B^2*C*H^2*R*z*3i - C*F*T*Y*z*1i - C*F*M*Y*z*1i + B^2*F*H^2*T*z - C*F*Y*z*3i - B^2*F*H^2*z + B^4*H^4*R^2*z - C^2*T*z + F^2*M*Y^2*z - 2*C^2*z + F^2*Y^2*z - 2*B^2*C*F*H^2*R*Y - B*C*F*H*M*Y*1i + B^3*F*H^3*M*R*Y - B^2*F^2*H^2*Y*1i + B^2*C^2*H^2*R*2i - B^4*C*H^4*R^2 + B*F^2*H*M*Y^2 - B^4*F*H^4*R*1i - B^2*C*F*H^2 + C^2*F*Y*2i - C*F^2*Y^2 + C^3;

% Use the coeffs function to obtain coefficients
coefficients = coeffs(polynomial, z);

% Display the coefficients
%disp('Coefficients:');
%disp(coefficients);
latex(coefficients)

