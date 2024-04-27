clc
clear
figure()
NU=linspace(40e-6,50e-6,5);
for m=1:5
    %nu=50e-6;%kinematic viscosity \nu=\mu/\rho
    nu=NU(m);
    g=9.8;
    sc=deg2rad(30);
    rint=2.5e-3;
    q0=3e-5/(2*pi);
    aint=30*rint;

    hs=((3*nu*q0/(g*sin(sc)))*(1/aint))^(1/3);
    us=(1/3)*(g*sin(sc)/nu)*(hs)^2;
    Ys=3*nu/(hs)^2;
    h2=aint*cos(sc);

    D=us;
    E=hs;
    F=cos(sc);
    G=g;
    H=1/h2;
    P=0;
    Q=6*nu/(hs)^3;
    
    %Y=-Ys;
    N1=2000;

    A=linspace(0,2000,N1);
    B=1 ;
    %x1=[];
    for j=1:N1
            Y=-Ys.*(1+(hs.^2).*(A(j))^2);
            SQRT=sqrt(Y^2-4*E*F*G*(B)^2*H^2-4i*D*E*Q*A(j)-4*E*F*G*(A(j))^2);
            %x1(j)=Y/2-real(SQRT)/2;
            x2(j)=Y/2+real(SQRT)/2;
    end
    %x1=Y/2-real(SQRT)/2;
    %x2=Y/2+real(SQRT)/2;

    strB{m}=[num2str(nu)];
    
    plot(A,x2./A)
    legend(strB)
    hold on
end
