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
wend = 10;

Nw = 100;

Omega = linspace(wint,wend,Nw);
I = eye(2);
Ze = zeros(2,2);
% without long-wave approx.

for i = 1:Nw
    omega = Omega(i);
    
    M211 = 0;
    M212 = 0;
    M221 = 0;
    M222 = 1/Re;
%     M222 = 0; % long-wave approx.
    
    M111 = 1i;
    M112 = 1i;
    M121 = 1i/Fr^2;
    M122 = 1i;
    
    M011 = -1i*omega;
    M012 = 0;
    M021 = -2*C2/Re;
    M022 = -1i*omega+C2/Re;
    
    M2 = [M211,M212;M221,M222];
    M1 = [M111,M112;M121,M122];
    M0 = [M011,M012;M021,M022];
    
    LHSnew = [M0,M1;Ze,I];
    RHSnew = [Ze,M2;-I,Ze];
    
    [V,Dnew] = eig(LHSnew,RHSnew);
    Klist(i,:) = -diag(Dnew);
    
end

    
% subplot(3,1,1)
% figure()
plot(Omega,-imag(Klist(:,1)),'-o')
hold on
plot(Omega,-imag(Klist(:,2)),'-o')
hold on
plot(Omega,-imag(Klist(:,3)),'-o')
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



