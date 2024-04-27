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

uw = 1.978613807;
% uw = 1.98;


[x,y] = ode45(@ostorder,[0 1000],[1.0]);

figure()
plot(x,y,'-o')
hold on

function dy = ostorder(x,y)
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

% uw = 1.978613807;
% uw = 1.98;
% uw = 1.9788;
uw = 2.03;

    dy = zeros(1,1); % x1 = y ; x2 = y'
    term4 = tan(sc1)+(tan(sc2)-tan(sc1))*(1-uw+uw*y(1))/(1-uw+uw*y(1)+gamma*y(1)^(5/2));
    dy(1) = (tan(sc)-term4)/(1-(F^2)*((uw-1)^2)/(y(1)^3));
    
%     term1 = (1/(2*y(1)))*(y(2)^2);
%     term2 = R*((y(1))^(3/2))/((F^2)*(uw-1));
%     term3 = (1-(F^2)*((uw-1)^2)/(y(1)^3))*y(2);
     
end