clc
clear

g = 9.8;
sc = deg2rad(4.6);

Re = 23*2/3;
nu = 4.89e-6;

hsinf = (3*nu^2*Re/(g*sin(sc)))^(1/3);
usinf = g*hsinf^2*sin(sc)/(3*nu);

Fr = usinf/((g*hsinf*cos(sc))^(1/2));

% LHS X = i\omega X

kint = 1e-6;
kend = 2;
Nk = 400;

kk1 = linspace(kint,kend,Nk);  

Omegalist1 = [];
Omegalist2 = [];
VlistH = [];
VlistU = [];

C2 = 3;

for i = 1:Nk
    k1 = kk1(i);
    
    LHS11 = 1i*k1;
    LHS12 = 1i*k1;
    LHS21 = 1i*k1/Fr-(2*C2)/Re;
    LHS22 = 1i*(1)*k1+(1/Re)*(k1^2)+C2/Re;
%     LHS22 = 1i*(1)*k1+3/Re;
    
    LHS = [LHS11,LHS12;LHS21,LHS22];
%     LHS = LHS./1i;
    [V,Dnew] = eig(LHS);
    Dnew = Dnew./1i;
    DiagD = diag(Dnew);
    Omegalist1(i,1) = DiagD(1,1);
    Omegalist2(i,1) = DiagD(2,1);
    
end
subplot(3,1,1)
plot(kk1',imag(Omegalist1)/(2*pi),'-o')
hold on
% plot(kk1',imag(Omegalist2),'-^')

subplot(3,1,2)

% plot(kk1',real(Omegalist1),'-o')
% hold on
plot(kk1,(real(Omegalist1)),'-o')
hold on

subplot(3,1,3)
plot(kk1,real(Omegalist1)./kk1','-^')
hold on
% figure()
% plot(kk1,'o')

%%
kai = deg2rad(60);
Cg = 7;

kr = kk1' - imag(Omegalist1).*sin(kai)./Cg;
ki = imag(Omegalist1).*cos(kai)./Cg;

% plot(kr,ki,'-o')
% hold on
plot(kr,real(Omegalist1)./kr,'o')

% if imag(DiagD(1,1))>imag(DiagD(2,1))
%         Omegalist(i,1) = DiagD(1,1);
%         VlistH(i,1) = V(1,1);
%         VlistU(i,1) = V(2,1);
%     elseif imag(DiagD(2,1))>imag(DiagD(1,1))
%         Omegalist(i,1) = DiagD(2,1);
%         VlistH(i,1) = V(1,2);
%         VlistU(i,1) = V(2,2);
%     else
%         print('Wrong!')
% end