clc
clear
nu=40e-6;%kinematic viscosity \nu=\mu/\rho
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
T=1.2;
%Y=-Ys;
%Y=-Ys.*(1+(hs.^2).*((A(j))^2+(H^2)*(B)^2));
ThreeDre=[];% store k2 series in 3Dplot

N1=100;
N2=100;
AA=linspace(1,500,N1);
BB=linspace(0,0,N2);
%A=100;
%B=1;
%figure()
for k=1:N2
    B=BB(k);
    x1=[];
    x2=[];
    x3=[];
    for j=1:N1
        A=AA(j);
        
        Y=-Ys.*(1+(1/3).*(hs.^2).*((A)^2+(H^2)*(B)^2));
            aaa=1;
            %bbb=-(2.*Y-A.*D.*3i).*1i;
            bbb=-(2.*Y-T.*A.*D.*3i).*1i;
            %ccc=(-(A.^2).*(D^2).*3i+(Y.^2).*1i+4.*A.*D.*Y-A.*D.*E.*Q+(A.^2).*E*F*G*1i+(B.^2).*E*F*G*(H^2)*1i).*1i;
            ccc=(-(1+2.*T).*(A.^2).*(D^2).*1i+(Y.^2).*1i+(3+T).*A.*D.*Y-A.*D.*E.*Q+(A.^2).*E*F*G*1i+(B.^2).*E*F*G*(H^2)*1i).*1i;
            ddd=-(A.^3).*(D^3).*T+A.*D*(Y.^2)-(1+T).*(A.^2).*(D^2).*Y.*1i+(A.^2).*(D^2)*E*Q*1i+(A.^3).*D*E*F*G+(A.^2).*E*F*G.*Y.*1i-A.*D*E*Q.*Y+(B.^2).*E*F*G*(H^2).*Y.*1i+(T-1).*A*(B^2)*D*E*F*G*(H^2);

       
        coeffs=[aaa,bbb,ccc,ddd];
        result=roots(coeffs);
        x1(j)=imag(result(1,1));
        x2(j)=imag(result(2,1));
        x3(j)=imag(result(3,1));
    end
    ThreeDre(k,:)=x3;
subplot(3,1,1)
plot(AA,x1)
hold on
subplot(3,1,2)
plot(AA,x2)
hold on
subplot(3,1,3)
plot(AA,x3)
hold on
xx3=x3';
end
% figure()
% surf(AA.*hs,BB.*H.*hs,ThreeDre)
% figure()
% subplot(3,1,1)
% surf(AA.*hs,BB.*H.*hs,ThreeDre)
% view(0,90);
% title('Top view')
% subplot(3,1,2)
% surf(AA.*hs,BB.*H.*hs,ThreeDre)
% view(-90,0);
% title('Left view (k_1 fixed)')
% subplot(3,1,3)
% surf(AA.*hs,BB.*H.*hs,ThreeDre)
% view(0,0);
% title('Top view (k_2 fixed)')
% colormap winter

% figure()
% [Graph,c]=contour(AA.*hs,BB.*H.*hs,ThreeDre,[0,0],'red');
% hold on
% c.LineWidth = 3;
% 
% [Graph,c]=contour(AA.*hs,BB.*H.*hs,ThreeDre,20);
% colormap winter
% c.LineWidth = 3;

%AAA=AA';
%disp(result);



