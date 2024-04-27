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

%clc
N1 = 401;
N2 = 401;

K1 = linspace(0,1,N1);
K2 = linspace(0,0.4,N2);
Ntry = 2001;
TTRY = linspace(1,2001,Ntry);
Nint = 5;
Nend = 71;
Ninterval = 1;
for zeta = Nint:Ninterval:Nend
%TRY = 41; % get the location of base flow
TRY = TTRY(zeta);
Hs = HS(TRY);
Us = US(TRY);
dHs = dHS(TRY);
dUs = dUS(TRY);
gammastar = gamma*(Hs^(3/2))/Us;
mus = tan(sc1)+(tan(sc2)-tan(sc1))/(1+gammastar);
Gampstar = (tan(sc2)-tan(sc1))*gammastar/(1+gammastar)^2;
Ain = A(TRY);
for i = 1:N2
    k2 = K2(i);
    for j = 1:N1
        k1 = K1(j);
        
        L11 = 1i*k1*Us+dUs+1*Us/Ain; % curvature term Us/A
        L12 = 1i*k1*Hs+dHs+1*Hs/Ain; % curvature term Hs/A
        L13 = 1i*k2*Hs;
        L21 = (1/F^2)*(1i*k1-(3/2)*(Gampstar/Hs));
        L22 = 1i*k1*Us+dUs+(1/R)*(Hs^(1/2))*(k1^2+k2^2)+(Gampstar/(F^2))/Us;
        L23 = 0;
        L31 = 1i*k2/(F^2);
        L32 = 0;
        L33 = 1i*k1*Us+(1/R)*(Hs^(1/2))*(k1^2+k2^2)+(1/F^2)*(mus/Us);
        
        LHS = [L11,L12,L13;
            L21,L22,L23;
            L31,L32,L33];
        RHS = 1i.*eye(3);
        
        [Vnew,Dnew] = eig(LHS,RHS);
        Eigenvaltemp = diag(Dnew);
        Eigenval(:,j) = Eigenvaltemp; % record all eigenvalue for each k2 round 
    end
Growthomega = imag(Eigenval);
Realomega = real(Eigenval);

% subplot(2,1,1)
% plot(K1,imag(Eigenval(1,:)),'-o')
% hold on
% plot(K1,imag(Eigenval(2,:)),'-o')
% hold on
% plot(K1,imag(Eigenval(3,:)),'-o')
% hold on

% plot(K1,max(Growthomega),'-','linewidth',2)
% hold on
% xlabel('k_1')
% ylabel('\omega_i')
% 
% subplot(2,1,2)
% plot(A,HS,'o')
% line([0,2000],[1,1],'linestyle','--','linewidth',2,'color','red');
% hold on
% plot(A(TRY),Hs,'hexagram','Markersize',15,'MarkerFaceColor','red')
% xlabel('A')
% ylabel('H_s')

OneG(i,:) = Growthomega(1,:);
TwoG(i,:) = Growthomega(2,:);
ThreeG(i,:) = Growthomega(3,:);
TotalG(i,:) = max(Growthomega);

end

% Find the max value and location in contour
% location1->K2 ; location2->K1
[OneGtemp,OneGtemplocation] = max(OneG);
[OneGmax,OneGlocation2] = max(OneGtemp);
OneGlocation1 = OneGtemplocation(OneGlocation2);

[TwoGtemp,TwoGtemplocation] = max(TwoG);
[TwoGmax,TwoGlocation2] = max(TwoGtemp);
TwoGlocation1 = TwoGtemplocation(TwoGlocation2);

[ThreeGtemp,ThreeGtemplocation] = max(ThreeG);
[ThreeGmax,ThreeGlocation2] = max(ThreeGtemp);
ThreeGlocation1 = ThreeGtemplocation(ThreeGlocation2);

[TotalGtemp,TotalGtemplocation] = max(TotalG);
[TotalGmax,TotalGlocation2] = max(TotalGtemp);
TotalGlocation1 = TotalGtemplocation(TotalGlocation2);

% maxlevel whole surf
% subplot(1,2,1)
% [Graph,c]=contour(K1,K2,TotalG,[0,max(TotalG)]);
% hold on
% [Graphbb,cbb]=contour(K1,K2,TotalG,'LevelList',[0.01,0.01],'LineColor','red','LineWidth',1,'ShowText','on');
% 
% c.LineWidth = 1;
% xlabel('k_1')
% ylabel('k_2')
% legendstr=['Res=',num2str(Res),'   ','Fr=',num2str(Fr)];
% legend(legendstr)
% colorbar('northoutside')
% colormap winter
% 
% subplot(1,2,2)
% plot(A,HS,'o')
% line([0,2000],[1,1],'linestyle','--','linewidth',2,'color','red');
% hold on
% plot(A(TRY),Hs,'hexagram','Markersize',15,'MarkerFaceColor','red')
% xlabel('A')
% ylabel('H_s')

% max-levellist comparison
% subplot(1,2,1)
% [Graphbb,cbb]=contour(K1,K2,TotalG,'linecolor','r','LevelList',[0.01,0.01],'LineWidth',1,'ShowText','on');
% hold on
% c.LineWidth = 1;
% xlabel('k_1')
% ylabel('k_2')
% legendstr=['Res=',num2str(Res),'   ','Fr=',num2str(Fr)];
% title(legendstr)
% colormap winter

% subplot(1,2,2)
% plot(A,HS,'o')
% line([0,2000],[1,1],'linestyle','--','linewidth',2,'color','red');
% hold on
% plot(A(TRY),Hs,'hexagram','Markersize',15,'MarkerFaceColor','r')
% xlabel('A')
% ylabel('H_s')

% Dif-level solutions (1G,2G,3G,TotalG) comparison    

% subplot(3,2,1)
% [Graphbbb,cbbb]=contour(K1,K2,OneG,[0,max(OneG)]);
% hold on
% plot(K1(OneGlocation2),K2(OneGlocation1),'hexagram','Markersize',15,'MarkerFaceColor','red')
% c.LineWidth = 1;
% xlabel('k_1')
% ylabel('k_2')
% colorbar('northoutside')
% colormap winter
% 
% subplot(3,2,2)
% [Graphb4,cb4]=contour(K1,K2,TwoG,[0,max(TwoG)]);
% hold on
% plot(K1(TwoGlocation2),K2(TwoGlocation1),'hexagram','Markersize',15,'MarkerFaceColor','red')
% c.LineWidth = 1;
% xlabel('k_1')
% ylabel('k_2')
% colorbar('northoutside')
% colormap winter
% 
% subplot(3,2,3)
% [Graphb5,cb5]=contour(K1,K2,ThreeG,[0,max(ThreeG)]);
% hold on
% plot(K1(ThreeGlocation2),K2(ThreeGlocation1),'hexagram','Markersize',15,'MarkerFaceColor','red')
% c.LineWidth = 1;
% xlabel('k_1')
% ylabel('k_2')
% colorbar('northoutside')
% colormap winter
% 
% subplot(3,2,4)
% [Graphb6,cb6]=contour(K1,K2,TotalG,[0,max(TotalG)]);
% hold on
% plot(K1(TotalGlocation2),K2(TotalGlocation1),'hexagram','Markersize',15,'MarkerFaceColor','red')
% c.LineWidth = 1;
% xlabel('k_1')
% ylabel('k_2')
% colorbar('northoutside')
% colormap winter
% 
% subplot(3,2,[5,6])
% plot(A,HS,'o')
% line([0,2000],[1,1],'linestyle','--','linewidth',2,'color','red');
% hold on
% plot(A(TRY),Hs,'hexagram','Markersize',15,'MarkerFaceColor','red')
% xlabel('A')
% ylabel('H_s')
% 
% legendstr=['Res=',num2str(Res),'   ','Fr=',num2str(Fr)];
% title(legendstr)

% Dif-level solutions (1G,2G,3G,TotalG) maxpoint movement
subplot(1,2,1)
plot(K1(TotalGlocation2),K2(TotalGlocation1),'hexagram','Markersize',10,'MarkerFaceColor','red')
hold on
c.LineWidth = 1;
xlabel('k_1')
ylabel('k_2')
colorbar('northoutside')
colormap winter
subplot(1,2,2)
plot(A,HS,'o')
line([0,2000],[1,1],'linestyle','--','linewidth',2,'color','red');
hold on
plot(A(TRY),Hs,'hexagram','Markersize',10,'MarkerFaceColor','red')
xlabel('A')
ylabel('H_s')

legendstr=['Res=',num2str(Res),'   ','Fr=',num2str(Fr)];
title(legendstr)


end


% r,g,b,m colorbar list
% useless tools
% quiver(A(TRY),HS(TRY),0,1,'MaxHeadSize',20,'Linestyle','-','Linewidth',0.5,'color','b')
% hold on
% text(A(TRY),HS(TRY),'\leftarrow')
% MapA = rescale(A,0,1);
% MapHS = rescale(HS,0,1);




























