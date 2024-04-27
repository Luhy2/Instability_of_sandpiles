clc
clear
%nu=40e-6;%kinematic viscosity \nu=\mu/\rho
NUM=2000;
NU=linspace(4e-6,50e-6,NUM);

N1=4000;
N2=5;
AA=linspace(0,1000,N1);
BB=linspace(1,2,N2);
%A=100;
%B=1;
figure()
for k=1:N2
    B=BB(k);
    
    pd=[]; % pd and z are used to find zero points under dif conditions
    z=1; % Count num
    
    x3=[];
    for m=1:NUM
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
        %Y=-Ys.*(1+(hs.^2).*((A(j))^2+(H^2)*(B)^2));
        for j=1:N1
            A=AA(j);

            Y=-Ys.*(1+(hs.^2).*((A)^2+(H^2)*(B)^2));

            aaa=1;
            bbb=-(2.*Y-A.*D.*3i).*1i;
            ccc=(-(A.^2).*(D^2).*3i+(Y.^2).*1i+4.*A.*D.*Y-A.*D.*E.*Q+(A.^2).*E*F*G*1i+(B.^2).*E*F*G*(H^2)*1i).*1i;
            ddd=-(A.^3).*(D^3)+A.*D*(Y.^2)-(A.^2).*(D^2).*Y.*2i+(A.^2).*(D^2)*E*Q*1i+(A.^3).*D*E*F*G+(A.^2).*E*F*G.*Y.*1i-A.*D*E*Q.*Y+(B.^2).*E*F*G*(H^2).*Y.*1i;

            coeffs=[aaa,bbb,ccc,ddd];
            result=roots(coeffs);
            %x1(j)=imag(result(1,1));
            %x2(j)=imag(result(2,1));
            x3(j)=imag(result(3,1));
        end
        for n= 1:N1-1

            %if x2(n)==0
                %pd(1,z)=(A(n)+A(n+1))/2;
                %pd(2,z)=nu;
                %z=z+1;
            %elseif x2(n)*x2(n+1)<0
            if x3(n)*x3(n+1)<0
                  pd(1,z)=(AA(n)+AA(n+1))/2;
                  pd(2,z)=nu;
                  z=z+1;
            end
        end
        
        
    end
    
    % plot for phase diagram of k - nu
    Re=(hs*us)./pd(2,:);
    plot(Re,pd(1,:))
    
    strB{k}=['k_2=',num2str(B)];
    legend(strB)
    
    xlabel('Re')
    ylabel('k_1')
    
    hold on
    
    
end
%AAA=AA';
%disp(result);



