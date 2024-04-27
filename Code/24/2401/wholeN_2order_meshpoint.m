clc
clear

% Origin input parameters
Fr=1.02;
Res=8.45;
Gam=0.05877;

% Nspace=101;%400
Nspace=201;%300
% Nspace default 200
Aint=1100;
Aend=1210; %nondimension length by rint -> hsref
A=linspace(Aint,Aend,Nspace);
h=(Aend-Aint)/(Nspace-1);%each interval length

sc1=deg2rad(20.9);
sc2=deg2rad(32.76);
sc=deg2rad(29); %倾角
charL=0.825e-3;
beta=0.136;
g=9.8;
%C=5/4; %profile shape coeff.
    nu=(2/9)*(charL/beta)*sqrt(g)*(sin(sc)/sqrt(cos(sc)))*((tan(sc2)-tan(sc))/(tan(sc)-tan(sc1)));
    gammastar=(tan(sc2)-tan(sc))/(tan(sc)-tan(sc1));
    gamma=(beta/charL)*(g*cos(sc))^(1/2);

% Downstream info.
    usref=(Fr*Res*nu*(g*cos(sc))^(1/2))^(1/2);
    hsref=(Res/Fr)*nu/(g*cos(sc))^(1/2);
    Rint=hsref;
    aend=Aend*Rint;
    a=A.*Rint; %real a
    %Q/2pi here determined by init Fr and Res from 'infty' downstream
    q0=aend*cos(sc)*hsref*usref;
    hs=(q0*gammastar./(gamma.*a)).^(2/5);
    us=(gamma/gammastar).*hs.^(3/2);%gammastar是文章里的gamma

%use downstream as ref value of h and u
rint=hsref;
uint=usref;

%New nondimen Numbers
F0=uint/(g*rint*cos(sc))^(1/2);
R0=(uint*(rint)^(1/2))/nu;
% us,hs maybe changed when not steady uniform base flow
Us=us./uint;
Hs=hs./rint;
    
    G=gammastar;
    eta=0.01*tan(sc1)+(tan(sc2)-tan(sc1))/(1+G);
    Gamp=(G./((1+G).^2)).*(tan(sc2)-tan(sc1));
%     eta=0;
 % Mark coeffs.
 F=F0;
 R=R0;
 C=5/4;
 %B default 2.6 ~ K2 0.01
 B=2.6;
 X2=1/cos(sc);
 
%  plot(A,Us,'o')
%  hold on
 
 %%
 
 
% Examine Hs-A basic flow
%  plot(A,Hs./A,'-o')
%  hold on
L1=Us./(2.*h);
L2=Hs./(2.*h);
L3=Us./A;
L4=Hs./A;
L5=1i.*B.*X2.*Hs./A;
L6=-Us./(2.*h);
L7=-Hs./(2.*h);

L8=1./(2.*(F.^2).*h).*ones(1,Nspace);

L9=(C.*Us./(2.*h))-(1/R).*(Hs.^(1/2))./(h.^2);
L10=-(3/2).*(Gamp./(F.^2))./Hs;
L11=2.*Us./A+(Gamp./(F.^2))./Us+(1/R).*(Hs.^(1/2)).*(2./(h.^2))+(1/R).*(Hs.^(1/2)).*((B.*X2).^2)./(A.^2);
L12=1i.*B.*X2.*Us./A;

L13=-1./(2.*(F.^2).*h).*ones(1,Nspace);

L14=-(C.*Us./(2.*h))-(1/R).*(Hs.^(1/2))./(h.^2);
L15=Us./(2.*h)-(1/R).*(Hs.^(1/2))./(h.^2);
L16=1i.*(1/(F.^2)).*B.*X2./A;
L17=(eta./F.^2)./Us+(1/R).*(Hs.^(1/2)).*(2./(h.^2))+(1/R).*(Hs.^(1/2)).*((B.*X2).^2)./(A.^2);
L18=-Us./(2.*h)-(1/R).*(Hs.^(1/2))./(h.^2);

%倒数第三行的h系数
M1=Us(1,Nspace)./(2.*h);
M2=Hs(1,Nspace)./(2.*h);
M3=(-4.*Us(1,Nspace))./(2.*h);
M4=(-4.*Hs(1,Nspace))./(2.*h);
M5=((3.*Us(1,Nspace))./(2.*h))+Us(1,Nspace)./A(1,Nspace);
M6=((3.*Hs(1,Nspace))./(2.*h))+Hs(1,Nspace)./A(1,Nspace);
M7=1i.*B.*X2.*Hs(1,Nspace)./A(1,Nspace);

K1=3./(2.*h);
K2=-2./h;
K3=1./(2.*h);

% K1 = 1;
% K2 = 0;
% K3 = 0;


% init LHS and RHS
MATSIZE=(Nspace).*3;
LHS=zeros(MATSIZE);
rHS=zeros(MATSIZE);

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
                LHS(i,1)=1;
            end
            if lkg==2
                LHS(i,2)=1;
            end
            if lkg==3
                LHS(i,3)=1;
            end
        
        elseif lg==Nspace %right B.C.s
            if lkg==1
%                 LHS(i,3*lg-2)=K1;
%                 LHS(i,3*(lg-1)-2)=K2;
%                 LHS(i,3*(lg-2)-2)=K3;
                  
                LHS(i,3*lg)=M7;
                LHS(i,3*lg-1)=M6;
                LHS(i,3*lg-2)=M5;
                LHS(i,3*(lg-1)-1)=M4;
                LHS(i,3*(lg-1)-2)=M3;
                LHS(i,3*(lg-2)-1)=M2;
                LHS(i,3*(lg-2)-2)=M1;
            end
            if lkg==2
                LHS(i,3*lg-1)=K1;
                LHS(i,3*(lg-1)-1)=K2;
                LHS(i,3*(lg-2)-1)=K3;
            end
            if lkg==3
                LHS(i,3*lg)=K1;
                LHS(i,3*(lg-1))=K2;
                LHS(i,3*(lg-2))=K3;
            end
        else % interior points
            if lkg==1
                LHS(i,3*(lg+1)-1)=L2(1,lg);
                LHS(i,3*(lg+1)-2)=L1(1,lg);
                LHS(i,3*lg)=L5(1,lg);
                LHS(i,3*lg-1)=L4(1,lg);
                LHS(i,3*lg-2)=L3(1,lg);
                LHS(i,3*(lg-1)-1)=L7(1,lg);
                LHS(i,3*(lg-1)-2)=L6(1,lg);
            end
            if lkg==2
                LHS(i,3*(lg+1)-1)=L9(1,lg);
                LHS(i,3*(lg+1)-2)=L8(1,lg);
                LHS(i,3*lg)=L12(1,lg);
                LHS(i,3*lg-1)=L11(1,lg);
                LHS(i,3*lg-2)=L10(1,lg);
                LHS(i,3*(lg-1)-1)=L14(1,lg);
                LHS(i,3*(lg-1)-2)=L13(1,lg);
            end
            if lkg==3
                LHS(i,3*(lg+1))=L15(1,lg);
                LHS(i,3*lg)=L17(1,lg);
                LHS(i,3*lg-2)=L16(1,lg);
                LHS(i,3*(lg-1))=L18(1,lg);
            end
        end
    end
    lkg=lkg+1;
end
        
for i=1:MATSIZE
    for j=1:MATSIZE
        if i==j
            if i>3
                if i<3*Nspace-1
                    rHS(i,j)=1;
                end
            end
        end
    end
end

RHS=1i.*rHS;
%%
[V,D]=eig(LHS,RHS);
% [V,D]=eigs(LHS,RHS,3*Nspace-5,'smallestabs');
% [V,D]=eigs(LHS,RHS,50,'smallestabs');
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

% hold off;

%%
% figure()
TRY=472;
for i=1:MATSIZE/3
    AlphaH(i,1)=V((i-1)*3+1,TRY);
    AlphaU(i,1)=V((i-1)*3+2,TRY);
    AlphaV(i,1)=V((i-1)*3+3,TRY);
end
% plot(Ahalf,AlphaH,'-o')
% figure()
% plot(Ahalf,AlphaU)
% figure()

figure()
subplot(3,1,1)
plot(A(1:end),real(AlphaH),'-o')
hold on
plot(A(1:end),imag(AlphaH))
hold on

subplot(3,1,2)
plot(A(1:end),real(AlphaU),'-o')
hold on
plot(A(1:end),imag(AlphaU))
hold on

subplot(3,1,3)
plot(A(1:end),real(AlphaV),'-o')
hold on
plot(A(1:end),imag(AlphaV))
hold on

%%
AbsU=abs(AlphaH);
figure()

semilogy(A,AbsU,'-o')

%%


% verification of LHS mat

% for i=1:MATSIZE
%     if lkg==4
%         lkg=1;
%         lg=lg+1;
%     end
%     for j=1:MATSIZE
%         if lg==1 %left B.C.s
%             if lkg==1
%                 LHS(i,1)=1;
%             end
%             if lkg==2
%                 LHS(i,2)=1;
%             end
%             if lkg==3
%                 LHS(i,3)=1;
%             end
%         
%         elseif lg==Nspace %right B.C.s
%             if lkg==1
%                 LHS(i,3*lg-2)=1000;
%                 LHS(i,3*(lg-1)-2)=2000;
%                 LHS(i,3*(lg-2)-2)=3000;
%             end
%             if lkg==2
%                 LHS(i,3*lg-1)=1000;
%                 LHS(i,3*(lg-1)-1)=2000;
%                 LHS(i,3*(lg-2)-1)=3000;
%             end
%             if lkg==3
%                 LHS(i,3*lg)=1000;
%                 LHS(i,3*(lg-1))=2000;
%                 LHS(i,3*(lg-2))=3000;
%             end
%         else % interior points
%             if lkg==1
%                 LHS(i,3*(lg+1)-1)=2;
%                 LHS(i,3*(lg+1)-2)=1;
%                 LHS(i,3*lg)=5;
%                 LHS(i,3*lg-1)=4;
%                 LHS(i,3*lg-2)=3;
%                 LHS(i,3*(lg-1)-1)=7;
%                 LHS(i,3*(lg-1)-2)=6;
%             end
%             if lkg==2
%                 LHS(i,3*(lg+1)-1)=9;
%                 LHS(i,3*(lg+1)-2)=8;
%                 LHS(i,3*lg)=12;
%                 LHS(i,3*lg-1)=11;
%                 LHS(i,3*lg-2)=10;
%                 LHS(i,3*(lg-1)-1)=14;
%                 LHS(i,3*(lg-1)-2)=13;
%             end
%             if lkg==3
%                 LHS(i,3*(lg+1))=15;
%                 LHS(i,3*lg)=17;
%                 LHS(i,3*lg-2)=16;
%                 LHS(i,3*(lg-1))=18;
%             end
%         end
%     end
%     lkg=lkg+1;
% end









 
 