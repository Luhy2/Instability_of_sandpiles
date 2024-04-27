clc
clear

    rint=2.5e-3;
    Aint=1;
    Aend=400;
    Num=500;
    A=linspace(Aint,Aend,Num);
    AAA=linspace(1,400,500);
    a=AAA.*rint;
    
    
    sc1=deg2rad(20.9);
    sc2=deg2rad(32.76);
    sc=deg2rad(29); %倾角
    charL=0.825e-3;
    beta=0.136;
    g=9.8;

    nu=(2/9)*(charL/beta)*sqrt(g)*(sin(sc)/sqrt(cos(sc)))*((tan(sc2)-tan(sc))/(tan(sc)-tan(sc1)));

    q0=3.5e-6;%Q/2pi
    uint=2*q0/rint^2;
    
    gammastar=(tan(sc2)-tan(sc))/(tan(sc)-tan(sc1));
    gamma=(beta/charL)*(g*cos(sc))^(1/2);
    lambdaB=gamma/q0;
    
    %暂时用参数化表达hs与us
    hs=(q0*gammastar./(gamma.*a)).^(2/5);
    us=(gamma/gammastar).*hs.^(3/2);%gammastar是文章里的gamma
    
    
    % 便捷计算代号
    C=5/4; %Velocity profile shape coeff.
    
    %从后向前计算
    yprime=@(x,y)((g.*cos(sc).*(tan(sc)-tan(sc1)-(tan(sc2)-tan(sc1))./(1+lambdaB.*x.*y.^(5/2))).*(y.^3).*(x.^3)+C*(q0^2).*y+2*q0*nu.*y.^(5/2))./(-C*(q0^2).*x+g*cos(sc).*(x.^3).*(y.^3)-2*q0*nu.*x.*y.^(3/2)));
    
    xint=Aint*rint;
    xend=Aend*rint;
    
    xspan=[xend xint];
    hsint=(q0*gammastar./(gamma.*xend)).^(2/5);
    
    [x,y]=ode45(yprime,xspan,hsint);
    v=q0./(x.*y); %us
    
    %从前向后计算
%     Yprime=@(X,Y)((g.*cos(sc).*(tan(sc)-tan(sc1)-(tan(sc2)-tan(sc1))./(1+lambdaB.*X.*Y.^(5/2))).*(Y.^3).*(X.^3)+C*(q0^2).*Y+2*q0*nu.*Y.^(5/2))./(-C*(q0^2).*X+g*cos(sc).*(X.^3).*(Y.^3)-2*q0*nu.*X.*Y.^(3/2)));
%     
%     xint=Aint*rint;
%     xend=Aend*rint;
%     
%     Xspan=[xint xend];
%     Hsint=q0./(rint.*uint);
%     
%     [X,Y]=ode45(Yprime,Xspan,Hsint);
%     V=q0./(X.*Y); %us
    
    
    
    
    
    subplot(2,1,1)
    plot(A,hs./rint,'-o')
    hold on
    plot(x./rint,y./rint)
    hold on
%     plot(X./rint,Y./rint)
%     hold on
    
    subplot(2,1,2)
    plot(A,us./uint,'-o')
    hold on
    plot(x./rint,v./uint)
    hold on
%     plot(X./rint,V./uint)
%     hold on
    
    
    
    
    