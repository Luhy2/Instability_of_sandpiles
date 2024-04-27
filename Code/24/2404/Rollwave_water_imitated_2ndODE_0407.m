clc
clear

g=9.8;
sc=deg2rad(4.6);
rint=2.5e-3;

%q0=3e-5/(2*pi);
q0 = 3.534e-3/(2*pi);

aint = 80*rint;
aend = 2000*rint;

c0=q0;

nu=4.89e-6;

C=c0;
c1=27/35;
c2=3;

EXT=c2*nu*C;
hint=(C*c2*nu/(g*sin(sc)*aend))^(1/3);

N = 450;
% N = 650;

Nr = N+1;

a = linspace(aint,aend,Nr);
hs = (C*c2*nu./(g*sin(sc).*a)).^(1/3);
us = C./(a.*hs);
% test base flow
% subplot(2,1,1)
% plot(a,hs,'o')
% subplot(2,1,2)
% plot(a(1:end-1),diff(hs),'o')
% hold on
% plot(a,us,'o')

usinf = us(1,end);
hsinf = hs(1,end);

Us = us./usinf;
Hs = hs./hsinf;
A = a./hsinf;
% test nond base flow
% figure()
% subplot(2,1,1)
% plot(A,Hs,'o')
% subplot(2,1,2)
% plot(A(1:end-1),diff(Hs),'o')

R = usinf*hsinf/nu;
F = usinf/(g*hsinf*cos(sc))^(1/2);

%------------------------------------------------


uw = 2.273;
% uw = 1.98;


[x,y] = ode45(@secorder,[80 1000],[1,0]);

subplot(2,1,1)
plot(x,y(:,1),'-o')
hold on
subplot(2,1,2)
% plot(x,y(:,2),'-o')
plot(y(:,1),y(:,2),'-o')
hold on



function dy = secorder(x,y)
g=9.8;
sc=deg2rad(4.6);
rint=2.5e-3;

%q0=3e-5/(2*pi);
q0 = 3.534e-3/(2*pi);

aint = 80*rint;
aend = 2000*rint;

c0=q0;

nu=4.89e-6;

C=c0;
c1=27/35;
c2=3;

EXT=c2*nu*C;
hint=(C*c2*nu/(g*sin(sc)*aend))^(1/3);

N = 450;
% N = 650;

Nr = N+1;

a = linspace(aint,aend,Nr);
hs = (C*c2*nu./(g*sin(sc).*a)).^(1/3);
us = C./(a.*hs);

usinf = us(1,end);
hsinf = hs(1,end);

Us = us./usinf;
Hs = hs./hsinf;
A = a./hsinf;

R = usinf*hsinf/nu;
F = usinf/(g*hsinf*cos(sc))^(1/2);
%---------------------------------------------
% uw = 1.978613807;
% uw = 1.98;
% uw = 1.9788;
uw = 2.2730;

    dy = zeros(2,1); % x1 = y ; x2 = y'
    dy(1) = y(2);
    term1 = (2/(y(1)))*(y(2)^2);
    term2 = R*((y(1))^(2))/((F^2)*(uw-1));
    term3 = (1-(F^2)*((uw-1)^2)/(y(1)^3))*y(2);
    term4 = ((F^2)/R)*(3/(y(1)^2))*(uw+(1-uw)/y(1));
    
    dy(2) = term1+term2*(term3-tan(sc)+term4);
    
end




