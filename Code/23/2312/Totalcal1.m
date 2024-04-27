clc
clear

% % 参数

sc1=deg2rad(20.9);
sc2=deg2rad(32.76);
sc=deg2rad(29); %倾角
charL=0.825e-3;
beta=0.136;
g=9.8;

nu=(2/9)*(charL/beta)*sqrt(g)*(sin(sc)/sqrt(cos(sc)))*((tan(sc2)-tan(sc))/(tan(sc)-tan(sc1)));

rint=2.5e-3;

aint=1*rint;
aend=50*rint;

Gam=0.05877;
q0=3.5e-7;%Q/2pi
uint=2*q0/rint^2;

%gammastar=(tan(sc2)-tan(sc))/(tan(sc)-tan(sc1));
%gamma=(beta/charL)*(g*cos(sc))^(1/2);
Res=uint*(rint)^(1/2)/nu;
Fr=uint/(g*cos(sc)*rint)^(1/2);


tint=aint/rint;
tend=aend/rint;
amesh=linspace(tint,tend,100);
solinit=bvpinit(amesh,@guess);

opts = bvpset('RelTol',1e-17,'Stats','on');%指定相对误差容限
sol=bvp5c(@aODE,@bcfun,solinit,opts);

%sol=bvp5c(@aODE,@bcfun,solinit);
%hs=(q0*gammastar./(gamma.*sol.x)).^(2/5);
%us=(gamma/gammastar).*hs.^(3/2);

%figure()
subplot(3,2,1)
plot(sol.x,real(sol.y(1,:)),'-o')
hold on
%legend('\omega=0','\omega=0.5','\omega=1','\omega=1.5','\omega=2','location','Best')
legend('k_2=0','k_2=10','k_2=50','k_2=100','k_2=200','location','Best')
xlabel('A')
yl=ylabel('\alpha_H');
%title('k_2=0.5 case')
title('\omega=0.5 case')
%set(yl,'Rotation',0)

subplot(3,2,2)
plot(sol.x,real(sol.y(2,:)),'-o')
hold on
legend('k_2=0','k_2=10','k_2=50','k_2=100','k_2=200','location','Best')
xlabel('A')
yl=ylabel('\alpha_U');
title('\omega=0.5 case')

subplot(3,2,3)
plot(sol.x,real(sol.y(3,:)),'-o')
hold on
legend('k_2=0','k_2=10','k_2=50','k_2=100','k_2=200','location','Best')
xlabel('A')
yl=ylabel('\alpha_V');
title('\omega=0.5 case')
resui=real(sol.y(3,:))./(sol.x);

 subplot(3,2,4)
 plot(sol.x,real(sol.y(3,:))./(sol.x),'-o')
 hold on
 legend('k_2=0','k_2=10','k_2=50','k_2=100','k_2=200','location','Best')
 xlabel('A')
 yl=ylabel('\alpha_V/A');
 title('\omega=0.5 case')
 
subplot(3,2,5)
plot(sol.x,real(sol.y(4,:)),'-o')
hold on
legend('k_2=0','k_2=10','k_2=50','k_2=100','k_2=200','location','Best')
xlabel('A')
yl=ylabel('d\alpha_U/dA');
title('\omega=0.5 case')

subplot(3,2,6)
plot(sol.x,real(sol.y(5,:)),'-o')
hold on
legend('k_2=0','k_2=10','k_2=50','k_2=100','k_2=200','location','Best')
xlabel('A')
yl=ylabel('d\alpha_V/dA');
title('\omega=0.5 case')

figure()
subplot(3,2,1)
plot(sol.x,imag(sol.y(1,:)),'-o')
hold on
subplot(3,2,2)
plot(sol.x,imag(sol.y(2,:)),'-o')
hold on
subplot(3,2,3)
plot(sol.x,imag(sol.y(3,:)),'-o')
hold on
subplot(3,2,5)
plot(sol.x,imag(sol.y(4,:)),'-o')
hold on
subplot(3,2,6)
plot(sol.x,imag(sol.y(5,:)),'-o')
hold on




