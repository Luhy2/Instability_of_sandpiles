clc
clear

g=9.8;
sc=deg2rad(20);
rint=2.5e-3;

%q0=3e-5/(2*pi);
q0 = 3e-3/(2*pi);

aint=80*rint;
aend=1500*rint;

c0=q0;

nu=50e-6;

C=c0;
c1=27/35;
c2=3;

EXT=c2*nu*C;
hint=(C*c2*nu/(g*sin(sc)*aend))^(1/3);

N = 400;

a = linspace(aint,aend,N);
hs = (C*c2*nu./(g*sin(sc).*a)).^(1/3);
us = C./(a.*hs);
% test base flow
% plot(a,hs,'o')
% hold on
% plot(a,us,'o')

usinf = us(1,end);
hsinf = hs(1,end);

R = usinf*hsinf/nu;
F = usinf/(g*hsinf*cos(sc))^(1/2);

%A = a./rint;
Aint = 800;
Aend = 1500;
A = linspace(Aint,Aend,N);
Ahalf = A(1:end-1)+diff(A)/2;

DeltaA = diff(A);

h = DeltaA(1,end);

C1 = 1; %1 2 14
C2 = 10;
C3 = 2400; % 10

% FA=1+C1*(1+tanh(C2*(-1+(2.*A./(Aint+Aend+C3)))));
FA=ones(1,N);
%%
% examine mesh and FA
% plot(A,'o')
% hold on
% plot(Ahalf,'o')
% figure()
plot(A,FA,'o')
hold on
%%

B = 0;

s1 = 1/(2*h);
s2 = 1/h;
s3 = (1/2)*1i*B./(Ahalf.*cos(sc));

s4 = -1/(2*h)-(1/(R*h^2)).*FA;
s5 = FA.*(-3/R)-1/(h*F^2);
s6 = (2/(R*h^2)+3/R).*FA;
s7 = FA.*(-3/R)+1/(h*F^2);
s8 = 1/(2*h)-(1/(R*h^2)).*FA;
s9 = (1/2)*(1/F^2)*1i*B./(A.*cos(sc));

K1=3./(2.*h);
K2=-2./h;
K3=1./(2.*h);

% K1=1;
% K2=-1;
% K3=0;

MATSIZE = 3*N;

LHS=zeros(MATSIZE);
rHS=eye(MATSIZE);

lg=1;%LHS counting switch times (current location)
lkg=1;%counting every 3 rows

for i=1:MATSIZE
    if lkg==4
        lkg=1;
        lg=lg+1;
    end
    for j=1:MATSIZE
        if lg==1 %left B.C.s
            if lkg==1
                LHS(i,1)=1/2;
                LHS(i,4)=1/2;
                rHS(i,1)=0;
            end
            if lkg==2
                LHS(i,2)=1;
                rHS(i,2)=0;
            end
            if lkg==3
                LHS(i,3)=1;
                rHS(i,3)=0;
            end
        elseif lg==N %right B.C.s
            if lkg==1
                LHS(i,3*lg)=s3(1,end);
                LHS(i,3*lg-1)=1/h;
                
                LHS(i,3*lg-2)=3/(2*h);%center
                
                LHS(i,3*lg-3)=s3(1,end);
                LHS(i,3*lg-4)=-1/h;
                LHS(i,3*lg-5)=-2/h;
                LHS(i,3*lg-8)=1/(2*h);
                
%                 LHS(i,3*lg-2)=1;
%                 rHS(i,3*lg-2)=0;
            end
            if lkg==2
%                 LHS(i,3*lg-1)=1;

%                 rHS(i,3*lg-1)=0;
%                 LHS(i,3*lg-1)=usinf.*K1-1i*usinf*0.01;
                LHS(i,3*lg-1)=usinf.*K1;
                LHS(i,3*(lg-1)-1)=usinf.*K2;
                LHS(i,3*(lg-2)-1)=usinf.*K3;
                
            end
            if lkg==3
%                 LHS(i,3*lg)=1;
%                 rHS(i,3*lg)=0;
%                 LHS(i,3*lg)=usinf.*K1-1i*usinf*0.01;
                LHS(i,3*lg)=usinf.*K1;
                LHS(i,3*(lg-1))=usinf.*K2;
                LHS(i,3*(lg-2))=usinf.*K3;
                
            end
        else % interior points
            if lkg==1
                LHS(i,3*lg+1)=s1;
                LHS(i,3*lg)=s3(1,lg-1);
                LHS(i,3*lg-1)=s2;
                
                LHS(i,3*lg-2)=0;%center
                
                LHS(i,3*lg-3)=s3(1,lg-1);
                LHS(i,3*lg-4)=-s2;
                LHS(i,3*lg-5)=-s1;
            end
            if lkg==2
                LHS(i,3*lg+2)=s8(1,lg);
                LHS(i,3*lg+1)=s7(1,lg);
                
                LHS(i,3*lg-1)=s6(1,lg);
                
                LHS(i,3*lg-2)=s5(1,lg);
                LHS(i,3*lg-4)=s4(1,lg);
            end
            if lkg==3
                LHS(i,3*lg+3)=s8(1,lg);
                LHS(i,3*lg+1)=s9(1,lg);
                
                LHS(i,3*lg)=s6(1,lg);
                
                LHS(i,3*lg-2)=s9(1,lg);
                LHS(i,3*lg-3)=s4(1,lg);
            end
        end
    end
    lkg=lkg+1;
end


RHS = 1i.*rHS;


[V,D]=eig(LHS,RHS);

Omegalist=diag(D);
%%
figure()
realOme=real(Omegalist);
imagOme=imag(Omegalist);
plot(realOme,imagOme,'o')
hold on
% legend('200','400','600')
for i = 1:1:length(Omegalist)
    text(realOme(i), imagOme(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end
% xlim([0,2])
% ylim([-1,0])
% hold off;

%%
% figure()
TRY=445;
for i=1:MATSIZE/3
    AlphaH(i,1)=V((i-1)*3+1,TRY);
    AlphaU(i,1)=V((i-1)*3+2,TRY);
    AlphaV(i,1)=V((i-1)*3+3,TRY);
end
% plot(Ahalf,AlphaH,'-o')
% figure()
% plot(Ahalf,AlphaU)
% figure()
AAhalf = [Ahalf(1,1)-h Ahalf];

figure()
subplot(3,1,1)
plot(AAhalf(1:end),real(AlphaH),'-o')
hold on
plot(AAhalf(1:end),imag(AlphaH))
hold on
% xlim([800,3900])


subplot(3,1,2)
plot(A(1:end),real(AlphaU),'-o')
hold on
plot(A(1:end),imag(AlphaU))
hold on
% xlim([800,3900])

subplot(3,1,3)
plot(A(1:end),real(AlphaV),'-o')
hold on
plot(A(1:end),imag(AlphaV))
hold on
% xlim([800,900])



