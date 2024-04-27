clc
clear

g = 9.8;
sc = deg2rad(4.6);

Re = 23*2/3;
nu = 4.89e-6;

C2 = 3;
hsinf = (C2*nu^2*Re/(g*sin(sc)))^(1/3);
usinf = g*hsinf^2*sin(sc)/(C2*nu);

Fr = usinf/((g*hsinf*cos(sc))^(1/2));

normw = 2*pi*hsinf/usinf;
% w = normw*f
% wint = normw*0; 
% wend = normw*7;

wint = 0;
wend = 1;

Nw = 10;

Omega = linspace(wint,wend,Nw);

% long-wave approx.

for i = 1:Nw
    omega = Omega(i);
    
    LHS11 = -1i*omega;
    LHS12 = 0;
    LHS21 = -2*C2/Re;
    LHS22 = -1i*omega+C2/Re;
    
    RHS11 = -1i;
    RHS12 = -1i;
    RHS21 = -1i/Fr;
    RHS22 = -1i;
    
    LHS = [LHS11,LHS12;LHS21,LHS22];
    RHS = [RHS11,RHS12;RHS21,RHS22];
    
    [V,Dnew] = eig(LHS,RHS);
    Klist(i,:) = diag(Dnew);
    
end

    
% subplot(3,1,1)
plot(Omega,imag(Klist(:,1)),'-o')
hold on
% plot(Omega,imag(Klist(:,2)),'-o')
% hold on

% figure()
% plot(real(Klist(:,2)),imag(Klist(:,2)),'o')
% hold on

% subplot(3,1,2)
% 
% 
% plot(kk1,(real(Omegalist1)),'-o')
% hold on
% 
% subplot(3,1,3)
% plot(kk1,real(Omegalist1)./kk1','-^')
% hold on



