clc
clear
figure()
Nkb=5;
B=linspace(1,2,Nkb);
for k=1:Nkb
    nu=40e-6;%kinematic viscosity \nu=\mu/\rho
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

    A=linspace(0,150,N1);
    %B=0;
    x1=[];
    for j=1:N1
            Y=-Ys.*(1+(hs.^2).*((A(j))^2+(H^2)*(B(k))^2));
            SQRT=(Y^2-4*E*F*G*(B(k))^2*H^2-4i*D*E*Q*A(j)-4*E*F*G*(A(j))^2)^(1/2);
            
            x1(j)=Y/2-real(SQRT)/2;
            x2(j)=Y/2+real(SQRT)/2;
    end
    %x1=Y/2-real(SQRT)/2;
    %x2=Y/2+real(SQRT)/2;

    %KK=((A).^2+(H.^2).*B.^2).^(1/2);
    strB{k}=['k_2=',num2str(B(k))];
    %subplot(2,1,1)
    title('different k_2')
    plot(A,x2)
    legend(strB);
    xlabel('k_1')
    ylabel('Im{\omega}')
    hold on
    
    %subplot(2,1,2)
    %title('different k_2')
    %plot(log(A),x2./A)
    %legend(strB);
    %xlabel('log{k_1}')
    %ylabel('Im{\omega}/k_1')
    %hold on
    
end
