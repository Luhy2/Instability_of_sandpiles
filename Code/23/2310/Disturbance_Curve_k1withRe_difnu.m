clc
clear
figure()
NUM=5;
NU=linspace(40e-6,50e-6,NUM);

NUMKB=1;
KB=linspace(1,2,NUMKB);
for kb=1:NUMKB
    pd=[]; % pd and z are used to find zero points under dif conditions
    z=1;
    for m=1:NUM
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

        A=linspace(0,300,N1);
        B=KB(kb) ;
        %x1=[];
        for j=1:N1
                Y=-Ys.*(1+(hs.^2).*((A(j))^2+(H^2)*(B)^2));
                SQRT=sqrt(Y^2-4*E*F*G*(B)^2*H^2-4i*D*E*Q*A(j)-4*E*F*G*(A(j))^2);
                %x1(j)=Y/2-real(SQRT)/2;
                x2(j)=Y/2+real(SQRT)/2;
        end
        %x1=Y/2-real(SQRT)/2;
        %x2=Y/2+real(SQRT)/2;

        % plot for IM\omega - k
        strB{m}=['\nu=',num2str(nu)];

        plot(A,x2)
        xlabel('k_1')
        ylabel('Im[\omega]')
        legend(strB)
        hold on



    end


end





