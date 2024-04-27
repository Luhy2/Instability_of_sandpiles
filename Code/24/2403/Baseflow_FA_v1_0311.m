clc
clear

% init para
g = 9.8;
sc = deg2rad(20);
rint = 2.5e-3;
Qstar = 3e-5; % init Q flow rate
qstar = Qstar/(2*pi);
uint = (2*qstar)/(rint^2);
q0 = qstar/(uint*(rint^2));



nu = 50e-6; % kinematic viscosity \nu=\mu/\rho

R0 = uint*rint/nu; % 2*R0 = Rej
F0 = sqrt((uint^2)/(g*rint*cos(sc)));

% FA = 1;

c1 = 2*27/35;
c2 = 3;

Aint = 40;
Aend = 10; % inverse cal

% call ode45 func
Hint = (((F0^2)/R0)*c2*q0/(Aint*tan(sc)))^(1/3);

    % y -> H ; x -> A
    yprime = @(x,y) (c1*(q0^2).*y+(1/F0^2)*tan(sc).*(y.^3).*(x.^3)-(1/R0)*c2*q0.*(x.^2))/((1/F0^2).*(y.^3).*(x.^3)-c1*(q0^2).*x);

% Aspan = [Aint Aend];
Nspace = 100;
Aspan = linspace(Aint,Aend,Nspace);

[x,y] = ode45(yprime,Aspan,Hint);

U = q0./(y.*x);

    % figure()
    plot(x,y,'o')
    hold on

%%
    % create ideal H(A) and U(A) distribution
    Nspace = 100;
    AA = linspace(Aint+140,Aend,Nspace);
    HA = (((F0^2)/R0)*c2*q0./(AA.*tan(sc))).^(1/3);
    UUA = q0./(HA.*AA);
    
    plot(AA,HA,'^')
    hold on
%     plot(AA,UUA,'^')
%     hold on
    
    % res plot between ideal and real H(A)
%     figure()
%     plot(AA,HA-y','-o')
%     hold on
%%
% use 'init' conidtion to solve ode45
clc

    YQprime = @(x,y)-((c1*(q0^2).*y+(1/F0^2)*tan(sc).*(y.^3).*(x.^3)-(1/R0)*c2*q0.*(x.^2))/((1/F0^2).*(y.^3).*(x.^3)-c1*(q0^2).*x));
Nspace = 100;
AsQ = flip(Aspan);

[X,Y] = ode45(YQprime,AsQ,y(end,1));
    plot(X,Y,'o')
    hold on
    
    
%%
clc
% Apply FA to \nu
Acint = 180;
Acend = 10;
Ncspace = 400;

Ac = linspace(Acint,Acend,Ncspace);
% plot(Ac,'o')
% hold on

Ac1 = Ac(1,1:end/2);
Ac2 = Ac(1,(end/2)+1:end);

AcM = Ac(1,end/2);
% plot(Ac1,'o')
% hold on 
% plot(Ac2,'o')
% hold on

% FA1 = ones(1,length(Ac2));
% FA2 = exp((Ac1-(AcM))./(AcM/4));
% FA = [FA2 FA1];

% plot(Ac,FA,'o')
% hold on

% NUc = nu.*FA; 

% Rc0 = uint*rint./NUc; % modified Rc0


%     Yprime = @(X,Y)(((c1*(q0^2).*Y+(1/F0^2)*tan(sc).*(Y.^3).*X.^3)-(1./Rc0).*c2*q0.*(X.^2))/((1/F0^2).*(Y.^3).*(X.^3)-c1*(q0^2).*X));
% Yprime = @(X,Y) piecefunc(X);
AcF = flip(Ac);

[X,Y] = ode45(@piecefunc,AcF,y(end,1));

UUc = q0./(Y.*X);

    % figure()
    plot(X,Y,'o')
    hold on
    xlabel('A')
    ylabel('H')
    legend('ideal','damp1.4','damp10')
    
%     plot(X,UUc,'o')
%     hold on
%     legend('ideal','damp1.4','damp10')
%     xlabel('A')
%     ylabel('U')
    
    function Yprime = piecefunc(X,Y)
    
        % init para
        g = 9.8;
        sc = deg2rad(20);
        rint = 2.5e-3;
        Qstar = 3e-5; 
        qstar = Qstar/(2*pi);
        uint = (2*qstar)/(rint^2);
        q0 = qstar/(uint*(rint^2));

        nu = 50e-6; 

        R0 = uint*rint/nu; 
        F0 = sqrt((uint^2)/(g*rint*cos(sc)));

        c1 = 2*27/35;
        c2 = 3;

        Acint = 180;
        Acend = 10;
        Ncspace = 100;
        Ac = linspace(Acint,Acend,Ncspace);

        AcM = Ac(1,end/2);
         
        if X <= AcM
            
            Yprime = -((c1*(q0^2).*Y+(1/F0^2)*tan(sc).*(Y.^3).*X.^3)-(1/R0).*c2*q0.*(X.^2))/((1/F0^2).*(Y.^3).*(X.^3)-c1*(q0^2).*X);
        else
            FA = exp((X-(AcM))./(AcM/10));
            NUc = nu.*FA;
            Rc0 = uint*rint./NUc;
            
%             Yprime = 0;
            Yprime = -((c1*(q0^2).*Y+(1/F0^2)*tan(sc).*(Y.^3).*X.^3)-(1/Rc0).*c2*q0.*(X.^2))/((1/F0^2).*(Y.^3).*(X.^3)-c1*(q0^2).*X);
        end
    end


