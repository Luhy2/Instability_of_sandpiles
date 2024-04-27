clc
clear

% Origin input parameters
Fr=1.02;
Res=8.45;
Gam=0.05877;

Nspace=1000;
% Nspace default 200
Aint=10;
Aend=410; %nondimension length by rint -> hsref
A=linspace(Aint,Aend,Nspace);
h=(Aend-Aint)/(Nspace-1);%each interval length

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
aend=Aend*Rint;
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
Ahalf=A(1:end-1)+diff(A)/2;
    ahalf=Ahalf.*Rint; %real a
    hshalf=(q0*gammastar./(gamma.*ahalf)).^(2/5);
    ushalf=(gamma/gammastar).*hshalf.^(3/2);
    Hshalf=hshalf./rint;
    Ushalf=ushalf./uint;
    %eta=(tan(sc2)-tan(sc1))/gammastar;
    
    [Numhalfrow,Numhalfcol]=size(Ahalf);%find half points column value
    
    G=gammastar;
    eta=tan(sc1)+(tan(sc2)-tan(sc1))/(1+G);
    Gamp=(G./((1+G).^2)).*(tan(sc2)-tan(sc1));
   
 % Mark coeffs.
 F=F0;
 R=R0;
 C=5/4;
 %B default 2.6 ~ K2 0.01
 B=2.6;
 
 % B=k2 maybe put it into cycle or set as a value
    L1=Ushalf./h+Ushalf./(2.*Ahalf);
    L2=-Ushalf./h+Ushalf./(2.*Ahalf);
    L3=Hshalf./(2.*Ahalf);
    L4=Hshalf./(2.*Ahalf);
    L5=1i*(B/cos(sc)).*Hshalf./(2.*Ahalf);
    L6=1i*(B/cos(sc)).*Hshalf./(2.*Ahalf);
    L7=Hshalf./2;
    L8=Hshalf./2;
    L9=1/((F^2)*h)-(3/2)*(Gamp/(F^2))./(2.*Hshalf);
    L10=-1/((F^2)*h)-(3/2)*(Gamp/(F^2))./(2.*Hshalf);
    L11=Ushalf./Ahalf+(Gamp/(F^2))./(2.*Ushalf)+(1/R).*(Hshalf.^(1/2)).*((B/cos(sc))^2)./(2.*(Ahalf.^2));
    L12=Ushalf./Ahalf+(Gamp/(F^2))./(2.*Ushalf)+(1/R).*(Hshalf.^(1/2)).*((B/cos(sc))^2)./(2.*(Ahalf.^2));
    L13=1i*(B/cos(sc)).*(Ushalf./(2.*Ahalf));
    L14=1i*(B/cos(sc)).*(Ushalf./(2.*Ahalf));
%     L13=ones(1,Numhalfcol).*0;
%     L14=ones(1,Numhalfcol).*0;
    L15=C.*Ushalf./2-(1/R).*(Hshalf.^(1/2))./h;
    L16=C.*Ushalf./2+(1/R).*(Hshalf.^(1/2))./h;
    L17=1i*(1/F^2)*(B/cos(sc))./(2.*Ahalf);
    L18=1i*(1/F^2)*(B/cos(sc))./(2.*Ahalf);
    L19=(eta/(F^2))./(2.*Ushalf)+(1/R).*(Hshalf.^(1/2)).*((B/cos(sc))^2)./(2.*(Ahalf.^2));
    L20=(eta/(F^2))./(2.*Ushalf)+(1/R).*(Hshalf.^(1/2)).*((B/cos(sc))^2)./(2.*(Ahalf.^2));
    L21=Ushalf./2-(1/R).*(Hshalf.^(1/2))./h;
    L22=Ushalf./2+(1/R).*(Hshalf.^(1/2))./h;
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

lg=0;%LHS counting switch times
lkg=1;%counting every 5 rows
for i=1:MATSIZE
    if lkg==6
        lkg=1;
        lg=lg+1;
    end
    for j=1:MATSIZE
        if lg==0 % left B.Cs
            if lkg==1
                LHS(i,1+5*lg)=L8(1,lg+1);
                LHS(i,3+5*lg)=L1(1,lg+1);
                LHS(i,4+5*lg)=L3(1,lg+1);
                LHS(i,5+5*lg)=L5(1,lg+1);
                LHS(i,6+5*lg)=L7(1,lg+1);
            end
            if lkg==2
                LHS(i,1+5*lg)=L16(1,lg+1);
                LHS(i,3+5*lg)=L9(1,lg+1);
                LHS(i,4+5*lg)=L11(1,lg+1);
                LHS(i,5+5*lg)=L13(1,lg+1);
                LHS(i,6+5*lg)=L15(1,lg+1);
            end
            if lkg==3
                LHS(i,2+5*lg)=L22(1,lg+1);
                LHS(i,3+5*lg)=L17(1,lg+1);
                LHS(i,5+5*lg)=L19(1,lg+1);
                LHS(i,7+5*lg)=L21(1,lg+1);
            end
            if lkg==4
                LHS(i,1+5*lg)=L26(1,lg+1);
                LHS(i,4+5*lg)=L23(1,lg+1);
                LHS(i,6+5*lg)=L25(1,lg+1);
            end
            if lkg==5
                LHS(i,2+5*lg)=L30(1,lg+1);
                LHS(i,5+5*lg)=L27(1,lg+1);
                LHS(i,7+5*lg)=L29(1,lg+1);
            end
        
        elseif lg==(MATSIZE/5)-1 % right B.Cs
            if lkg==1
                LHS(i,-2+5*lg)=L2(1,lg+1);
                LHS(i,-1+5*lg)=L4(1,lg+1);
                LHS(i,0+5*lg)=L6(1,lg+1);
                LHS(i,1+5*lg)=L8(1,lg+1);
                LHS(i,3+5*lg)=L1(1,lg+1);
                LHS(i,4+5*lg)=L3(1,lg+1);
                LHS(i,5+5*lg)=L5(1,lg+1);
            end
            if lkg==2
                LHS(i,-2+5*lg)=L10(1,lg+1);
                LHS(i,-1+5*lg)=L12(1,lg+1);
                LHS(i,0+5*lg)=L14(1,lg+1);
                LHS(i,1+5*lg)=L16(1,lg+1);
                LHS(i,3+5*lg)=L9(1,lg+1);
                LHS(i,4+5*lg)=L11(1,lg+1);
                LHS(i,5+5*lg)=L13(1,lg+1);
            end
            if lkg==3
                LHS(i,-2+5*lg)=L18(1,lg+1);
                LHS(i,0+5*lg)=L20(1,lg+1);
                LHS(i,2+5*lg)=L22(1,lg+1);
                LHS(i,3+5*lg)=L17(1,lg+1);
                LHS(i,5+5*lg)=L19(1,lg+1);
            end
            if lkg==4
                LHS(i,-1+5*lg)=L24(1,lg+1);
                LHS(i,1+5*lg)=L26(1,lg+1);
                LHS(i,4+5*lg)=L23(1,lg+1);
            end
            if lkg==5
                LHS(i,0+5*lg)=L28(1,lg+1);
                LHS(i,2+5*lg)=L30(1,lg+1);
                LHS(i,5+5*lg)=L27(1,lg+1);
            end
        else % interior points
            if lkg==1
                LHS(i,-2+5*lg)=L2(1,lg+1);
                LHS(i,-1+5*lg)=L4(1,lg+1);
                LHS(i,0+5*lg)=L6(1,lg+1);
                LHS(i,1+5*lg)=L8(1,lg+1);
                LHS(i,3+5*lg)=L1(1,lg+1);
                LHS(i,4+5*lg)=L3(1,lg+1);
                LHS(i,5+5*lg)=L5(1,lg+1);
                LHS(i,6+5*lg)=L7(1,lg+1);
            end
            if lkg==2
                LHS(i,-2+5*lg)=L10(1,lg+1);
                LHS(i,-1+5*lg)=L12(1,lg+1);
                LHS(i,0+5*lg)=L14(1,lg+1);
                LHS(i,1+5*lg)=L16(1,lg+1);
                LHS(i,3+5*lg)=L9(1,lg+1);
                LHS(i,4+5*lg)=L11(1,lg+1);
                LHS(i,5+5*lg)=L13(1,lg+1);
                LHS(i,6+5*lg)=L15(1,lg+1);
            end
            if lkg==3
                LHS(i,-2+5*lg)=L18(1,lg+1);
                LHS(i,0+5*lg)=L20(1,lg+1);
                LHS(i,2+5*lg)=L22(1,lg+1);
                LHS(i,3+5*lg)=L17(1,lg+1);
                LHS(i,5+5*lg)=L19(1,lg+1);
                LHS(i,7+5*lg)=L21(1,lg+1);
            end
            if lkg==4
                LHS(i,-1+5*lg)=L24(1,lg+1);
                LHS(i,1+5*lg)=L26(1,lg+1);
                LHS(i,4+5*lg)=L23(1,lg+1);
                LHS(i,6+5*lg)=L25(1,lg+1);
            end
            if lkg==5
                LHS(i,0+5*lg)=L28(1,lg+1);
                LHS(i,2+5*lg)=L30(1,lg+1);
                LHS(i,5+5*lg)=L27(1,lg+1);
                LHS(i,7+5*lg)=L29(1,lg+1);
            end
        end
    end
    lkg=lkg+1;
end

[V,D]=eig(LHS,RHS);

Omegalist=diag(D);

% default para settings
Omeganew=Omegalist(1:end-2*Numhalfcol);
% Vnew=V(:,1:end-2*Numhalfcol);

% Omeganew=Omegalist(Numhalfcol:end-2*Numhalfcol);
Vnew=V(:,1:end-2*Numhalfcol);
%%

realOme=real(Omeganew);
imagOme=imag(Omeganew);
plot(realOme,imagOme,'o')
hold on
% legend('200','400','800')
% for i = 1:1:length(Omeganew)
%     text(realOme(i), imagOme(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
% end
% 
% hold off;
%%
figure()
%clc
TRY=230;
for i=1:MATSIZE/5
    AlphaH(i,1)=Vnew((i-1)*5+1+2,TRY);
    AlphaU(i,1)=Vnew((i-1)*5+2+2,TRY);
    AlphaV(i,1)=Vnew((i-1)*5+3+2,TRY);
    PU(i,1)=Vnew((i-1)*5+1,TRY);
    PV(i,1)=Vnew((i-1)*5+2,TRY);
end
% plot(Ahalf,AlphaH,'-o')
% figure()
% plot(Ahalf,AlphaU)
% figure()

subplot(5,1,1)
plot(A(2:end),real(AlphaH),'-o')
hold on
plot(A(2:end),imag(AlphaH))
hold on

subplot(5,1,2)
plot(A(2:end),real(AlphaU),'-o')
hold on
plot(A(2:end),imag(AlphaU))
hold on

subplot(5,1,3)
plot(A(2:end),real(AlphaV),'-o')
hold on
plot(A(2:end),imag(AlphaV))
hold on


subplot(5,1,4)
plot(A(1:end-1),real(PU),'-o')
hold on
plot(A(1:end-1),imag(PU))
hold on

subplot(5,1,5)
plot(A(1:end-1),real(PV),'-o')
hold on
plot(A(1:end-1),imag(PV))
hold on

% Back up Check for LHS

% lg=0;%LHS counting switch times
% lkg=1;%counting every 5 rows
% for i=1:MATSIZE
%     if lkg==6
%         lkg=1;
%         lg=lg+1;
%     end
%     for j=1:MATSIZE
%         if lg==0 % left B.Cs
%             if lkg==1
%                 LHS(i,1+5*lg)=8;
%                 LHS(i,3+5*lg)=1;
%                 LHS(i,4+5*lg)=3;
%                 LHS(i,5+5*lg)=5;
%                 LHS(i,6+5*lg)=7;
%             end
%             if lkg==2
%                 LHS(i,1+5*lg)=16;
%                 LHS(i,3+5*lg)=9;
%                 LHS(i,4+5*lg)=11;
%                 LHS(i,5+5*lg)=13;
%                 LHS(i,6+5*lg)=15;
%             end
%             if lkg==3
%                 LHS(i,2+5*lg)=22;
%                 LHS(i,3+5*lg)=17;
%                 LHS(i,5+5*lg)=19;
%                 LHS(i,7+5*lg)=21;
%             end
%             if lkg==4
%                 LHS(i,1+5*lg)=26;
%                 LHS(i,4+5*lg)=23;
%                 LHS(i,6+5*lg)=25;
%             end
%             if lkg==5
%                 LHS(i,2+5*lg)=30;
%                 LHS(i,5+5*lg)=27;
%                 LHS(i,7+5*lg)=29;
%             end
%         
%         elseif lg==(MATSIZE/5)-1 % right B.Cs
%             if lkg==1
%                 LHS(i,-2+5*lg)=2;
%                 LHS(i,-1+5*lg)=4;
%                 LHS(i,0+5*lg)=6;
%                 LHS(i,1+5*lg)=8;
%                 LHS(i,3+5*lg)=1;
%                 LHS(i,4+5*lg)=3;
%                 LHS(i,5+5*lg)=5;
%             end
%             if lkg==2
%                 LHS(i,-2+5*lg)=10;
%                 LHS(i,-1+5*lg)=12;
%                 LHS(i,0+5*lg)=14;
%                 LHS(i,1+5*lg)=16;
%                 LHS(i,3+5*lg)=9;
%                 LHS(i,4+5*lg)=11;
%                 LHS(i,5+5*lg)=13;
%             end
%             if lkg==3
%                 LHS(i,-2+5*lg)=18;
%                 LHS(i,0+5*lg)=20;
%                 LHS(i,2+5*lg)=22;
%                 LHS(i,3+5*lg)=17;
%                 LHS(i,5+5*lg)=19;
%             end
%             if lkg==4
%                 LHS(i,-1+5*lg)=24;
%                 LHS(i,1+5*lg)=26;
%                 LHS(i,4+5*lg)=23;
%             end
%             if lkg==5
%                 LHS(i,0+5*lg)=28;
%                 LHS(i,2+5*lg)=30;
%                 LHS(i,5+5*lg)=27;
%             end
%         else % interior points
%             if lkg==1
%                 LHS(i,-2+5*lg)=2;
%                 LHS(i,-1+5*lg)=4;
%                 LHS(i,0+5*lg)=6;
%                 LHS(i,1+5*lg)=8;
%                 LHS(i,3+5*lg)=1;
%                 LHS(i,4+5*lg)=3;
%                 LHS(i,5+5*lg)=5;
%                 LHS(i,6+5*lg)=7;
%             end
%             if lkg==2
%                 LHS(i,-2+5*lg)=10;
%                 LHS(i,-1+5*lg)=12;
%                 LHS(i,0+5*lg)=14;
%                 LHS(i,1+5*lg)=16;
%                 LHS(i,3+5*lg)=9;
%                 LHS(i,4+5*lg)=11;
%                 LHS(i,5+5*lg)=13;
%                 LHS(i,6+5*lg)=15;
%             end
%             if lkg==3
%                 LHS(i,-2+5*lg)=18;
%                 LHS(i,0+5*lg)=20;
%                 LHS(i,2+5*lg)=22;
%                 LHS(i,3+5*lg)=17;
%                 LHS(i,5+5*lg)=19;
%                 LHS(i,7+5*lg)=21;
%             end
%             if lkg==4
%                 LHS(i,-1+5*lg)=24;
%                 LHS(i,1+5*lg)=26;
%                 LHS(i,4+5*lg)=23;
%                 LHS(i,6+5*lg)=25;
%             end
%             if lkg==5
%                 LHS(i,0+5*lg)=28;
%                 LHS(i,2+5*lg)=30;
%                 LHS(i,5+5*lg)=27;
%                 LHS(i,7+5*lg)=29;
%             end
%         end
%     end
%     lkg=lkg+1;
% end

    
    