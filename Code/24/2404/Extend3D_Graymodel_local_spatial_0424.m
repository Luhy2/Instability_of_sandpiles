clc
clear
% Origin input parameters
Fr=1.02;
Res=8.45;
% Gam=0.05877;
Ainf=2000; %infty location
sc1=deg2rad(20.9);
sc2=deg2rad(32.76);
sc=deg2rad(29); %倾角
charL=0.825e-3;
beta=0.136;
g=9.8;

    nu=(2/9)*(charL/beta)*sqrt(g)*(sin(sc)/sqrt(cos(sc)))*((tan(sc2)-tan(sc))/(tan(sc)-tan(sc1)));
    gamma=(tan(sc2)-tan(sc))/(tan(sc)-tan(sc1));
% Downstream info.
    u0=(Fr*Res*nu*(g*cos(sc))^(1/2))^(1/2);
    h0=(Res/Fr)*nu/(g*cos(sc))^(1/2);
    lammb = u0/(h0^(3/2));
% Incoming flow char.
Nbase = 2001;
ainf=Ainf*h0; %real a

aint = 80*h0;

a = linspace(aint,ainf,Nbase);
A = a./h0;
% h = diff(A);
deltaA = (A(end)-A(1))/(Nbase-1); %iso-spacing grids

q0=ainf*h0*u0;

%New nondimen Numbers
F=u0/(g*h0*cos(sc))^(1/2);
R=(u0*(h0)^(1/2))/nu;

hs = (q0./(lammb.*a)).^(2/5);
us = q0./(a.*hs);

HS=hs/h0;
US=us/u0;
% get the derivation of HS and US
for p = 1:Nbase
    if p == 1
        dHS(p) = (-3*HS(p)+4*HS(p+1)-HS(p+2))/(2*deltaA);
        dUS(p) = (-3*US(p)+4*US(p+1)-US(p+2))/(2*deltaA);
    elseif p == Nbase
        dHS(p) = (3*HS(p)-4*HS(p-1)+HS(p-2))/(2*deltaA);
        dUS(p) = (3*US(p)-4*US(p-1)+US(p-2))/(2*deltaA);
    else
        dHS(p) = (HS(p+1)-HS(p-1))/(2*deltaA);
        dUS(p) = (US(p+1)-US(p-1))/(2*deltaA);
    end
end

% verify base flow
% subplot(3,2,1)
% plot(A,HS,'o')
% line([0,2000],[1,1],'linestyle','--','linewidth',2,'color','red');
% subplot(3,2,2)
% plot(A,US,'o')
% line([0,2000],[1,1],'linestyle','--','linewidth',2,'color','red');
% subplot(3,2,3)
% plot(A,dHS,'o')
% line([0,2000],[0,0],'linestyle','--','linewidth',2,'color','red');
% subplot(3,2,4)
% plot(A,dUS,'o')
% line([0,2000],[0,0],'linestyle','--','linewidth',2,'color','red');
% subplot(3,2,6)
% plot(A,us./hs.^(3/2),'o')
% line([0,2000],[lammb,lammb],'linestyle','--','linewidth',2,'color','red');

Nbeta = 401;
Nomega = 401;

K2 = linspace(0,1,Nbeta);
Omega = linspace(0,4,Nomega);

Ntry = 2001;
TTRY = linspace(1,2001,Ntry);
Nint = 21;
Nend = 21;
Ninterval = 600;

for zeta = Nint:Ninterval:Nend
    TRY = TTRY(zeta);
    Hs = HS(TRY);
    Us = US(TRY);
    dHs = dHS(TRY);
    dUs = dUS(TRY);
    gammastar = gamma*(Hs^(3/2))/Us;
    mus = tan(sc1)+(tan(sc2)-tan(sc1))/(1+gammastar);
    Gampstar = (tan(sc2)-tan(sc1))*gammastar/(1+gammastar)^2;
    Ain = A(TRY);
    for i = 1:Nbeta
        k2 = K2(i);
        for j = 1:Nomega
            omega = Omega(j);
            
            M011 = -1i*omega+dUs+1*Us/Ain;
            M012 = dHs+1*Hs/Ain;
            M013 = 1i*k2*Hs;
            M021 = -(3/2)*(1/F^2)*Gampstar/Hs;
            M022 = -1i*omega+dUs+(1/R)*(Hs^(1/2))*(k2^2)+(Gampstar/(F^2))/Us;
            M023 = 0;
            M031 = 1i*k2/(F^2);
            M032 = 0;
            M033 = -1i*omega+(1/R)*(Hs^(1/2))*(k2^2)+(1/F^2)*(mus/Us);
            
            M111 = 1i*Us;
            M112 = 1i*Hs;
            M113 = 0;
            M121 = 1i/(F^2);
            M122 = 1i*Us;
            M123 = 0;
            M131 = 0;
            M132 = 0;
            M133 = 1i*Us;
            
            M211 = 0;
            M212 = 0;
            M213 = 0;
            M221 = 0;
            M222 = 1*(1/R)*(Hs^(1/2));
            M223 = 0;
            M231 = 0;
            M232 = 0;
            M233 = (1/R)*(Hs^(1/2));
            
            M0 = [M011,M012,M013;
                M021,M022,M023;
                M031,M032,M033];
            M1 = [M111,M112,M113;
                M121,M122,M123;
                M131,M132,M133];
            M2 = [M211,M212,M213;
                M221,M222,M223;
                M231,M232,M233];
            OneM = eye(3);
            Ze = zeros(3,3);
            
            LHS = [M0,M1;
                Ze,OneM];
            RHS = [Ze,-M2;
                OneM,Ze];
            
            [Vnew,Dnew] = eig(LHS,RHS);
            Eigenvaltemp = diag(Dnew);
            Eigenval(:,j) = Eigenvaltemp;
        end
        Growthalpha = -imag(Eigenval);
        Realalpha = real(Eigenval);
        
        OneG(i,:) = Growthalpha(1,:);
        TwoG(i,:) = Growthalpha(2,:);
        ThreeG(i,:) = Growthalpha(3,:);
        FourG(i,:) = Growthalpha(4,:);
        FiveG(i,:) = Growthalpha(5,:);
        SixG(i,:) = Growthalpha(6,:);
        for zk = 1:Nomega
            for row = 1:6
                if Growthalpha(row,zk)>0.1
                    Growthalpha(row,zk)=-100;
                    
                end
            end
        end
        TotalG(i,:) = max(Growthalpha);
    end

%     subplot(2,1,1)
%     plot(Omega,Growthalpha(1,:),'-o')
%     hold on
%     plot(Omega,Growthalpha(2,:),'-o')
%     hold on
%     plot(Omega,Growthalpha(3,:),'-o')
%     hold on
%     plot(Omega,Growthalpha(4,:),'-o')
%     hold on
%     plot(Omega,Growthalpha(5,:),'-o')
%     hold on
%     plot(Omega,Growthalpha(6,:),'-o')
%     hold on

%     plot(Omega,TotalG,'-o')
%     hold on
%     ylim([-0.02,0.01])
%     xlim([0,2])
%     subplot(2,1,2)
%     plot(A,HS,'o')
%     line([0,2000],[1,1],'linestyle','--','linewidth',2,'color','red');
%     hold on
%     plot(A(TRY),Hs,'hexagram','Markersize',15,'MarkerFaceColor','red')
%     xlabel('A')
%     ylabel('H_s')
    
%     plot(Omega,Omega./Realalpha(2,:),'-o')
%     hold on
%     plot(Omega,Omega./Realalpha(3,:),'-o')
%     hold on
%     ylim([1,3])
%     plot(Omega,max(Omega./Realalpha(2,:),Omega./Realalpha(3,:)),'-o')
%     hold on
    
end

% maxlevel whole surf
subplot(1,2,1)
[Graph,c] = contour(Omega,K2,TotalG,[0,max(TotalG)]);
hold on
[Graphbb,cbb]=contour(Omega,K2,TotalG,'LevelList',[0.01,0.01],'LineColor','red','LineWidth',1,'ShowText','on');

c.LineWidth = 1;
xlabel('\omega')
ylabel('k_2')
legendstr=['Res=',num2str(Res),'   ','Fr=',num2str(Fr)];
legend(legendstr)
colorbar('northoutside')
colormap winter
% 
subplot(1,2,2)
plot(A,HS,'o')
line([0,2000],[1,1],'linestyle','--','linewidth',2,'color','red');
hold on
plot(A(TRY),Hs,'hexagram','Markersize',15,'MarkerFaceColor','red')
xlabel('A')
ylabel('H_s')







