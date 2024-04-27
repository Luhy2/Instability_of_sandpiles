clc
clear

g = 9.8;
sc = deg2rad(29);
F = 1.02;

nu = 1.13e-3;
gamma = 0.517;
beta = 0.136;
CharL = 0.825e-3;
sc1 = deg2rad(20.9);
sc2 = deg2rad(32.76);

hsinf = CharL*gamma*F/beta;
usinf = F*(g*hsinf*cos(sc))^(1/2);

R = usinf*sqrt(hsinf)/nu;
Gamp = gamma*(tan(sc2)-tan(sc1))/(1+gamma)^2;

I = eye(2);
Ze = zeros(2,2);

wint = 0;
% wint = 0.995;
% wend = 0.262626;
% wend = 0.997;
% wend = 0.996316;
wend = 2;
Nw = 1000;
Omega = linspace(wint,wend,Nw);
% without long-wave approx.

for i = 1:Nw
    omega = Omega(i);
    
    M211 = 0;
    M212 = 0;
    M221 = 0;
    M222 = (F^2)/R;
%     M222 = 0; % long-wave approx.
    
    M111 = 1;
    M112 = 1;
    M121 = 1i;
    M122 = 1i*(F^2);
    
    M011 = -omega;
    M012 = 0;
    M021 = -(3/2)*Gamp;
    M022 = Gamp-1i*(F^2)*omega;
    
    M2 = [M211,M212;M221,M222];
    M1 = [M111,M112;M121,M122];
    M0 = [M011,M012;M021,M022];
    
    LHSnew = [M0,M1;Ze,I];
    RHSnew = [Ze,M2;-I,Ze];
    
    [V,Dnew] = eig(LHSnew,RHSnew);
    Klist(i,:) = -diag(Dnew);
    
end



% subplot(2,1,1)
% plot(Omega,-imag(Klist(:,1)),'-o')
% hold on
% plot(Omega,-imag(Klist(:,2)),'-o')
% hold on
% plot(Omega,-imag(Klist(:,3)),'-o')
% hold on
% 
% 
% line([0,2],[0,0],'linestyle','--');
% ylim([-0.03,0.02])
% 
% 
% subplot(2,1,2)
% plot(Omega,Omega'./real(Klist(:,1)),'-o')
% hold on
% plot(Omega,Omega'./real(Klist(:,2)),'-o')
% hold on
% plot(Omega,Omega'./real(Klist(:,3)),'-o')
% 
% ylim([1,3])

figure()

plot(real(Klist(:,1)),-imag(Klist(:,1)),'-o')
hold on
plot(real(Klist(:,2)),-imag(Klist(:,2)),'-o')
hold on
plot(real(Klist(:,3)),-imag(Klist(:,3)),'-o')
hold on


line([0,2],[0,0],'linestyle','--');
ylim([-0.03,0.02])



