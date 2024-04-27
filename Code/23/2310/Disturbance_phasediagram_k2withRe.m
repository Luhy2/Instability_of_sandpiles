clc
clear
figure()
NUM=500;
NU=linspace(10e-6,50e-6,NUM);

NUMKB=5;
KB=linspace(0,500,NUMKB);

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

        B=linspace(0,500,N1);
        A=KB(kb) ;
        %x1=[];
        for j=1:N1
                Y=-Ys.*(1+(hs.^2).*((A)^2+(H^2)*(B(j))^2));
                SQRT=sqrt(Y^2-4*E*F*G*(B(j))^2*H^2-4i*D*E*Q*(A.^2).^(1/2)-4*E*F*G*(A)^2);
                %x1(j)=Y/2-real(SQRT)/2;
                x2(j)=Y/2+real(SQRT)/2;
        end
        %x1=Y/2-real(SQRT)/2;
        %x2=Y/2+real(SQRT)/2;

        % plot for IM\omega - k
        %strB{m}=[num2str(nu)];

        %plot(A,x2)
        %legend(strB)
        %hold on

        % find zero point under dif \mu


        for n= 1:N1-1

            if x2(n)==0
                pd(1,z)=(B(n)+B(n+1))/2;
                pd(2,z)=nu;
                z=z+1;
            elseif x2(n)*x2(n+1)<0
                  pd(1,z)=(B(n)+B(n+1))/2;
                  pd(2,z)=nu;
                  z=z+1;
            end
        end
        
    end

    % plot for phase diagram of k - nu
    plot((hs*us)./pd(2,:),pd(1,:))
    
    strB{kb}=['k_1=',num2str(A)];
    legend(strB)
    
    xlabel('Re')
    ylabel('k_2')
    
    hold on
end





