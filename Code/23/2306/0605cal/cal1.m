clearvars -except QWEHG
clc

g=9.8;
sc=pi/6;
rint=0.01;%r=0.01m
q0=0.01;%q0=r*v0=0.01m*1m/s
c0=q0*rint*cos(sc);%c0=q0*r*cos\alpha

nu=1e-6;%kinematic viscosity \nu=\mu/\rho

C=c0;
EXT=2*nu*C;


yprime = @(x,y) (-(x.^3).*(y.^3).*g.*sin(sc)+y.*C.^2-EXT.*(x.^2))./((x.^3).*(y.^3).*g.*cos(sc)-x.*C.^2);

xint=rint*cos(sc);
xspan = [xint 0.4];
%initial h

%h0=0.004:0.001:0.006;
h0=0.005;

[x,y]=ode45(yprime,xspan,h0);

v = C./(x.*y); %u_a^*

y2=0.4.*(x-1)+10;
y3=9.8*sin(sc)./v;
%small quantities vis
%s1=y./x; %h/a
%s2=v./x; %u_a^*/a


figure()
subplot(2,1,1)
plot(x,y,'-o')
xlabel('a')
ylabel('h^*')

subplot(2,1,2)
plot(x,v,'-o')
xlabel('a')
ylabel('u_a^*')

%figure()
%plot(x,v,'-o')
%hold on
%plot(x,y2)
%xlabel('a')

%figure()
%plot(x,y3,'-o')


