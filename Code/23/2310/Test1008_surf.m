clc
clear

nu=50e-6;%kinematic viscosity \nu=\mu/\rho
g=9.8;
sc=deg2rad(30);
rint=2.5e-3;
q0=3e-5/(2*pi);
aint=30*rint;

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
Q=6*nu/(hs)^3;
Y=-Ys;

N1=200;
N2=200;
A=linspace(0,1000,N1);
B=linspace(0,1000,N2);
x1=[];
for j=1:N1
    for k=1:N2
        SQRT=sqrt(Y^2-4*E*F*G*(B(k))^2*H^2-4i*D*E*Q*A(j)-4*E*F*G*(A(j))^2);
        x1(j,k)=Y/2-real(SQRT)/2;
        x2(j,k)=Y/2+real(SQRT)/2;
        %if x2(j,k)<0
        %    x2(j,k)=nan;
        %end
    end
end
%x1=Y/2-real(SQRT)/2;
%x2=Y/2+real(SQRT)/2;

figure()
subplot(2,1,1)
surf(A,B,x1)
xlabel('k1')
ylabel('k2')
subplot(2,1,2)
surf(A,B,x2)
xlabel('k1')
ylabel('k2')
%figure()
%subplot(2,2,1)
%plot(A,x1)
%subplot(2,2,2)
%plot(A,x2)






