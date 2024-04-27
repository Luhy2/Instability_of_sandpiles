clc
clear

Aint = 200;
Am1 = 400;
Aend = 800;

N1 = 100;
N2 = 150; % dissipation mesh points -1

A1 = linspace(Aint,Am1,N1);
h1 = (Am1-Aint)/(N1-1);
Am2 = Am1+h1;

alpha = h1/(Am2^2); % increasing rate

A2(1) = Am2;

for i = 2:N2
    A2(i) = A2(i-1)*(1+alpha*A2(i-1));
end

A = [A1 A2];

Ahalf=A(1:end-1)+diff(A)/2;


C1 = 1;
C2 = 9;
C3 = 2000;

FA=1+C1*(1+tanh(C2*(-1+(2*A/(Aint+Am1+C3)))));

FAhalf=FA(1:end-1)+diff(FA)/2;
%% 
% Check mesh and FA
DeltaA = diff(A);

% plot(A,'o')
% hold on
% plot(Ahalf,'^')

% plot(DeltaA,'o')

plot(A,FA,'o')
hold on
plot(Ahalf,FAhalf,'^')
% legend('C3 = 500','C1 = 1300','C1 = 2000')
% plot(Ahalf,FAhalf,'^')
%%

h=h1;
Nspace = N1+N2;

% Origin input parameters
Fr=1.02;
Res=8.45;
Gam=0.05877;



sc1=deg2rad(20.9);
sc2=deg2rad(32.76);
sc=deg2rad(29); %倾角
charL=0.825e-3;
beta=0.136;
g=9.8;
C=5/4; %profile shape coeff.
    nu=(2/9)*(charL/beta)*sqrt(g)*(sin(sc)/sqrt(cos(sc)))*((tan(sc2)-tan(sc))/(tan(sc)-tan(sc1)));
    gammastar=(tan(sc2)-tan(sc))/(tan(sc)-tan(sc1));
    gamma=(beta/charL)*(g*cos(sc))^(1/2);

% Downstream info.
    usref=(Fr*Res*nu*(g*cos(sc))^(1/2))^(1/2);
    hsref=(Res/Fr)*nu/(g*cos(sc))^(1/2);
    
% Incoming flow char.
%Rint=2.5e-2;
Rint=hsref;
% aend=Aend*Rint;

aend=Am2*Rint;

a=A.*Rint; %real a
%q0=3.5e-7; %Q/2pi here determined by init Fr and Res from 'infty' downstream
q0=aend*cos(sc)*hsref*usref;
    hs=(q0*gammastar./(gamma.*a)).^(2/5);
    us=(gamma/gammastar).*hs.^(3/2);%gammastar是文章里的gamma

    
%uint=2*q0/(rint^2); %only ref velocity
%use downstream as ref value of h and u
rint=hsref;
uint=usref;

%New nondimen Numbers
F0=uint/(g*rint*cos(sc))^(1/2);
R0=(uint*(rint)^(1/2))/nu;
% us,hs maybe changed when not steady uniform base flow
Us=us./uint;
Hs=hs./rint;


% Half quantities definition
% Ahalf=A(1:end-1)+diff(A)/2; % see above block
    ahalf=Ahalf.*Rint; %real a
    hshalf=(q0*gammastar./(gamma.*ahalf)).^(2/5);
    ushalf=(gamma/gammastar).*hshalf.^(3/2);
    Hshalf=hshalf./rint;
    Ushalf=ushalf./uint;
    %eta=(tan(sc2)-tan(sc1))/gammastar;
    
    [Numhalfrow,Numhalfcol]=size(Ahalf);%find half points column value
    
    G=gammastar;
    eta=0.01.*(tan(sc1)+(tan(sc2)-tan(sc1))/(1+G));
    
    Gamp=(G./((1+G).^2)).*(tan(sc2)-tan(sc1));
    %eta=Gamp;
 % Mark coeffs.
 F=F0;
 R=R0;
 C=5/4;
 %B default 2.6 ~ K2 0.01
 B=1;
 
 % B=k2 maybe put it into cycle or set as a value
    L1=Ushalf./h+Ushalf./(2.*Ahalf);
    L2=-Ushalf./h+Ushalf./(2.*Ahalf);
    L3=Hshalf./(2.*Ahalf);
    L4=Hshalf./(2.*Ahalf);
    L5=1i*(B/cos(sc)).*Hshalf./(2.*Ahalf);
    L6=1i*(B/cos(sc)).*Hshalf./(2.*Ahalf);
    L7=Hshalf./2;
    L8=Hshalf./2;
    
    L9ori=1/((F^2)*h)-(3/2)*(Gamp/(F^2))./(2.*Hshalf);
    L9=1/((F^2)*h)+FAhalf.*(-(3/2)*(Gamp/(F^2))./(2.*Hshalf));
    
    L10ori=-1/((F^2)*h)-(3/2)*(Gamp/(F^2))./(2.*Hshalf);
    L10=-1/((F^2)*h)+FAhalf.*(-(3/2)*(Gamp/(F^2))./(2.*Hshalf)); 
   
    L11ori=Ushalf./Ahalf+(Gamp/(F^2))./(2.*Ushalf)+(1/R).*(Hshalf.^(1/2)).*((B/cos(sc))^2)./(2.*(Ahalf.^2));
    L11=Ushalf./Ahalf+FAhalf.*(Gamp/(F^2))./(2.*Ushalf)+FAhalf.*((1/R).*(Hshalf.^(1/2)).*((B/cos(sc))^2)./(2.*(Ahalf.^2)));
    
    L12ori=Ushalf./Ahalf+(Gamp/(F^2))./(2.*Ushalf)+(1/R).*(Hshalf.^(1/2)).*((B/cos(sc))^2)./(2.*(Ahalf.^2));
    L12=Ushalf./Ahalf+FAhalf.*(Gamp/(F^2))./(2.*Ushalf)+FAhalf.*((1/R).*(Hshalf.^(1/2)).*((B/cos(sc))^2)./(2.*(Ahalf.^2)));
    
    
    L13=1i*(B/cos(sc)).*(Ushalf./(2.*Ahalf));
    L14=1i*(B/cos(sc)).*(Ushalf./(2.*Ahalf));
%     L13=ones(1,Numhalfcol).*0;
%     L14=ones(1,Numhalfcol).*0;
    
    L15ori=C.*Ushalf./2-(1/R).*(Hshalf.^(1/2))./h;
    L15=C.*Ushalf./2+FAhalf.*(-(1/R).*(Hshalf.^(1/2))./h);
    
    L16ori=C.*Ushalf./2+(1/R).*(Hshalf.^(1/2))./h;
    L16=C.*Ushalf./2+FAhalf.*((1/R).*(Hshalf.^(1/2))./h);
    
    
    L17=1i*(1/F^2)*(B/cos(sc))./(2.*Ahalf);
    L18=1i*(1/F^2)*(B/cos(sc))./(2.*Ahalf);
    
    
    L19ori=(eta/(F^2))./(2.*Ushalf)+(1/R).*(Hshalf.^(1/2)).*((B/cos(sc))^2)./(2.*(Ahalf.^2));
    L19=FAhalf.*(eta/(F^2))./(2.*Ushalf)+FAhalf.*((1/R).*(Hshalf.^(1/2)).*((B/cos(sc))^2)./(2.*(Ahalf.^2)));
    
    L20ori=(eta/(F^2))./(2.*Ushalf)+(1/R).*(Hshalf.^(1/2)).*((B/cos(sc))^2)./(2.*(Ahalf.^2));
    L20=FAhalf.*(eta/(F^2))./(2.*Ushalf)+FAhalf.*((1/R).*(Hshalf.^(1/2)).*((B/cos(sc))^2)./(2.*(Ahalf.^2)));
    
    L21ori=Ushalf./2-(1/R).*(Hshalf.^(1/2))./h;
    L21=Ushalf./2+FAhalf.*(-(1/R).*(Hshalf.^(1/2))./h);
    
    L22ori=Ushalf./2+(1/R).*(Hshalf.^(1/2))./h;
    L22=Ushalf./2+FAhalf.*((1/R).*(Hshalf.^(1/2))./h);
    
    
    L23=ones(1,Numhalfcol).*(1/h);
    L24=ones(1,Numhalfcol).*(-1/h);
    L25=ones(1,Numhalfcol).*(-1/2);
    L26=ones(1,Numhalfcol).*(-1/2);
    L27=ones(1,Numhalfcol).*(1/h);
    L28=ones(1,Numhalfcol).*(-1/h);
    L29=ones(1,Numhalfcol).*(-1/2);
    L30=ones(1,Numhalfcol).*(-1/2);
    
% initialize LHS and RHS
MATSIZE=(Nspace-1)*5;
LHS=zeros(MATSIZE);
rHS=zeros(MATSIZE);

rg=1;%rHS counting
for i=1:MATSIZE
    for j=1:MATSIZE
        if rg==1
            if j==i+2
                if mod(i,5)<4
                    if mod(i,5)>0
                        rHS(i,j)=1;
                    end
                end
            end
        
        else
            if j==i+2
                if mod(i,5)<4
                    if mod(i,5)>0
                        rHS(i,j)=1;
                    end
                end
            end
            if j==i-3
                if mod(i,5)<4
                    if mod(i,5)>0
                        rHS(i,j)=1;
                    end
                end
            end
        end
    end
    if mod(i,5)==0
        rg=rg+1;
    end
end

RHS=(1i/2)*rHS;

lg=1;%LHS counting switch times
lkg=1;%counting every 5 rows
for i=1:MATSIZE
    if lkg==6
        lkg=1;
        lg=lg+1;
    end
    for j=1:MATSIZE+5
        
         % interior points
            if lkg==1
                LHS(i,-4+5*lg)=L2(1,lg);
                LHS(i,-3+5*lg)=L4(1,lg);
                LHS(i,-2+5*lg)=L6(1,lg);
                LHS(i,-1+5*lg)=L8(1,lg);
                LHS(i,1+5*lg)=L1(1,lg);
                LHS(i,2+5*lg)=L3(1,lg);
                LHS(i,3+5*lg)=L5(1,lg);
                LHS(i,4+5*lg)=L7(1,lg);
            end
            if lkg==2
                LHS(i,-4+5*lg)=L10(1,lg);
                LHS(i,-3+5*lg)=L12(1,lg);
                LHS(i,-2+5*lg)=L14(1,lg);
                LHS(i,-1+5*lg)=L16(1,lg);
                LHS(i,1+5*lg)=L9(1,lg);
                LHS(i,2+5*lg)=L11(1,lg);
                LHS(i,3+5*lg)=L13(1,lg);
                LHS(i,4+5*lg)=L15(1,lg);
            end
            if lkg==3
                LHS(i,-4+5*lg)=L18(1,lg);
                LHS(i,-2+5*lg)=L20(1,lg);
                LHS(i,0+5*lg)=L22(1,lg);
                LHS(i,1+5*lg)=L17(1,lg);
                LHS(i,3+5*lg)=L19(1,lg);
                LHS(i,5+5*lg)=L21(1,lg);
            end
            if lkg==4
                LHS(i,-3+5*lg)=L24(1,lg);
                LHS(i,-1+5*lg)=L26(1,lg);
                LHS(i,2+5*lg)=L23(1,lg);
                LHS(i,4+5*lg)=L25(1,lg);
            end
            if lkg==5
                LHS(i,-2+5*lg)=L28(1,lg);
                LHS(i,0+5*lg)=L30(1,lg);
                LHS(i,3+5*lg)=L27(1,lg);
                LHS(i,5+5*lg)=L29(1,lg);
            end
        end
    
    lkg=lkg+1;
end


% NEW expand matrix
RRHS=zeros(5*Nspace);
LLHS=zeros(5*Nspace);
LLHS(1,1)=1;
LLHS(2,2)=1;
LLHS(3,3)=1;

LLHS(4:end-2,:)=LHS;
RRHS(4:end-2,4:end-2)=RHS;
RRHS(4,1)=(1i/2);
RRHS(5,2)=(1i/2);
RRHS(6,3)=(1i/2);

% 右边界一阶导为零
% LLHS(5*Nspace,5*Nspace)=1;
% LLHS(5*Nspace-1,5*Nspace-1)=1;
% RRHS(5*Nspace,5*Nspace)=0;
% RRHS(5*Nspace-1,5*Nspace-1)=0;
% 右边界函数值为零
LLHS(5*Nspace,5*Nspace-2)=1;
LLHS(5*Nspace-1,5*Nspace-3)=1;
RRHS(5*Nspace,5*Nspace-2)=0;
RRHS(5*Nspace-1,5*Nspace-3)=0;


%%
[V,D]=eig(LLHS,RRHS);

% [V,D]=eigs(LLHS,RRHS,5*Nspace-5-2*Numhalfcol,'smallestabs');
% [V,D]=eigs(LLHS,RRHS,60,'smallestabs');

% [V,D]=eig(LHS(:,4:end-2),RHS);
% [V,D]=eigs(LHS(:,4:end-2),RHS,5*Nspace-5-2*Numhalfcol,'smallestabs');
Omegalist=diag(D);

% default para settings
% Omeganew=Omegalist(1:end-2*Numhalfcol-5);
Omeganew=Omegalist;
% Vnew=V(:,1:end-2*Numhalfcol);

% Omeganew=Omegalist(Numhalfcol:end-2*Numhalfcol);
% Vnew=V(:,1:end-2*Numhalfcol-5);
Vnew=V;

%%
% figure()
realOme=real(Omeganew);
imagOme=imag(Omeganew);
plot(realOme,imagOme,'o')
hold on
xlabel('\omega_R')
ylabel('\omega_I')
% legend('200','400','800')
for i = 1:1:length(Omeganew)
    text(realOme(i), imagOme(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end


% hold off;
%%
clc
%figure()
%clc
TRY=730;
for i=1:(MATSIZE/5)+1
    AlphaH(i,1)=Vnew((i-1)*5+1,TRY);
    AlphaU(i,1)=Vnew((i-1)*5+2,TRY);
    AlphaV(i,1)=Vnew((i-1)*5+3,TRY);
    PU(i,1)=Vnew((i-1)*5+4,TRY);
    PV(i,1)=Vnew((i-1)*5+5,TRY);
end

% AbsU=abs(AlphaU);
% figure()
% 
% % semilogy(A(1,1:200),real(AlphaU(1:200,1)),'-o')
% plot(A(1,1:200),real(AlphaU(1:200,1)),'-o')
% hold on
% plot(A(1,1:200),imag(AlphaU(1:200,1)),'-o')
% hold on

% for i=1:MATSIZE/5
%     AlphaH(i,1)=Vnew((i-1)*5+1+2,TRY);
%     AlphaU(i,1)=Vnew((i-1)*5+2+2,TRY);
%     AlphaV(i,1)=Vnew((i-1)*5+3+2,TRY);
%     PU(i,1)=Vnew((i-1)*5+1,TRY);
%     PV(i,1)=Vnew((i-1)*5+2,TRY);
% end

% plot(Ahalf,AlphaH,'-o')
% figure()
% plot(Ahalf,AlphaU)
% figure()
figure()

subplot(5,1,1)
plot(A(1:end),real(AlphaH),'-o')
hold on
plot(A(1:end),imag(AlphaH))
hold on
xlabel('A')
ylabel('Hprime')
% xlim([10,510])

subplot(5,1,2)
plot(A(1:end),real(AlphaU),'-o')
hold on
plot(A(1:end),imag(AlphaU))
hold on
xlabel('A')
ylabel('Uprime')
% xlim([10,510])

subplot(5,1,3)
plot(A(1:end),real(AlphaV),'-o')
hold on
plot(A(1:end),imag(AlphaV))
hold on
xlabel('A')
ylabel('Vprime')
% xlim([10,510])

subplot(5,1,4)
plot(A(1:end),real(PU),'-o')
hold on
plot(A(1:end),imag(PU))
hold on
xlabel('A')
ylabel('dUprime')
% xlim([10,510])

subplot(5,1,5)
plot(A(1:end),real(PV),'-o')
hold on
plot(A(1:end),imag(PV))
hold on
xlabel('A')
ylabel('dVprime')
% xlim([10,510])

    