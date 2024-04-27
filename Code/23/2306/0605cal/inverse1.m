clear
clc

g=9.8;
sc=deg2rad(20);
rint=2.5e-3;%r=0.01m
            %q0=1.7e-5/(pi*rint)*1.01;%q0=r*v0=0.01m*1m/s
            %c0=q0*rint*cos(sc);%c0=q0*r*cos\alpha
q0=3e-5/(2*pi);
%H0=4.8e-3;
%xint=rint*cos(sc);
xint=20*rint;
xend=6*rint;
aint=24*rint;
%v0=q0/(pi*rint^2);
%qprime=sqrt(v0^2+2*g*H0)/v0*q0;
%c0=0.5*sqrt(v0^2+2*g*H0)*rint^2;
c0=q0;
%h0=rint/2*cos(sc);

%h0=3*rint; %h0=4.8e-3;

nu=50e-6;%kinematic viscosity \nu=\mu/\rho

C=c0;
c1=27/35;
c2=3;

EXT=c2*nu*C;
hint=(C*c2*nu/(g*sin(sc)*aint))^(1/3);


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
subplot(3,1,1)
plot(x/rint,y/rint,'-o')
%ylim([0,3])
xlim([9,14])
xlabel('a/r_0')
ylabel('h/r_0')


subplot(3,1,2)
plot(x/rint,v,'-o')
xlim([9,14])
xlabel('a/r_0')
ylabel('u_a^*')

subplot(3,1,3)
plot(x,v,'-o')
xlim([9*rint,14*rint])
xlabel('a')
ylabel('u_a^*')
%figure()
%plot(x,v,'-o')
%hold on
%plot(x,y2)
%xlabel('a')

%figure()
%plot(x,y3,'-o')
