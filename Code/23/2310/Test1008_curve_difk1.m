clc
clear
figure()
Nkb=5;
A=linspace(0,5000,Nkb);
for k=1:Nkb
    nu=20e-6;%kinematic viscosity \nu=\mu/\rho
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

    N1=200;

    B=linspace(0,200,N1);
    %B=0;
    x1=[];
    for j=1:N1
            Y=-Ys.*(1+(hs.^2).*((A(k))^2+(H^2)*(B(j))^2));
            SQRT=(Y^2-4*E*F*G*(B(j))^2*H^2-4i*D*E*Q*A(k)-4*E*F*G*(A(k))^2)^(1/2);
            x1(j)=Y/2-real(SQRT)/2;
            x2(j)=Y/2+real(SQRT)/2;
    end
    %x1=Y/2-real(SQRT)/2;
    %x2=Y/2+real(SQRT)/2;


    strB{k}=['k_1=',num2str(A(k))];
    %subplot(2,1,1)
    %plot(B,x1)
    %legend(strB);
    %hold on
    
    %subplot(2,1,2)
    %plot(B,x2)
    %legend(strB);
    %hold on
    subplot(2,1,1)
    title('different k_1')
    plot(B,x2)
    legend(strB);
    xlabel('k_2')
    ylabel('Im{\omega_2}')
    hold on
    
    %subplot(2,1,2)
    %title('different k_1')
    %plot(B,x1)
    %legend(strB);
    %xlabel('k_2')
    %ylabel('Im{\omega_1}')
    %hold on
    
    
    %subplot(2,1,2)
    %title('different k_1')
    %plot(log(B),x2./B)
    %legend(strB);
    %xlabel('log{k_2}')
    %ylabel('Im{\omega}/k_2')
    %hold on
end
