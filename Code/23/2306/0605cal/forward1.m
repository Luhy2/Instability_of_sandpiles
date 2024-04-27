clear
clc

nu=50e-6;%kinematic viscosity \nu=\mu/\rho
g=9.8;
sc=deg2rad(20);
rint=2.5e-3;

q0=3e-5/(2*pi);
Q=q0*2*pi;

xbl=rint*0.3155*(Q/(nu*rint))^(1/3);

xint=1.426*xbl;
xend=10*rint;


c0=q0;




C=c0;
c1=27/35;
c2=3;

EXT=c2*nu*C;
%hint=(C*c2*nu/(g*sin(sc)*xint))^(1/3);
hint=2.308*rint*((nu*rint)/Q)^(1/3);
kk=hint/rint;


yprime = @(x,y) ((x.^3).*(y.^3).*g.*sin(sc)+2.*c1.*y.*C.^2-EXT.*(x.^2))./((x.^3).*(y.^3).*g.*cos(sc)-2.*c1.*x.*C.^2);


xspan = [xint xend];


[x,y]=ode45(yprime,xspan,hint);

v = C./(x.*y); %u_a^*

y2=0.4.*(x-1)+10;
y3=9.8*sin(sc)./v;
%small quantities vis
%s1=y./x; %h/a
%s2=v./x; %u_a^*/a

rx=x/rint;
ry=y/rint;
ResultMatrix=[rx ry];

figure()
subplot(2,1,1)
plot(x/rint,y/rint,'-o')
%ylim([0,3])
xlabel('a/r_0')
ylabel('h/r_0')


subplot(2,1,2)
plot(x/rint,v,'-o')
xlabel('a/r_0')
ylabel('u_a^*')

%figure()
%plot(x,v,'-o')
%hold on
%plot(x,y2)
%xlabel('a')

%figure()
%plot(x,y3,'-o')
