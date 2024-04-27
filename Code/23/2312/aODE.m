% function

function dy = aODEf(t,y)

    % 参数
    
    B=1;
    w=0.5;
    
    sc1=deg2rad(20.9);
    sc2=deg2rad(32.76);
    sc=deg2rad(29); %倾角
    charL=0.825e-3;
    beta=0.136;
    g=9.8;

    nu=(2/9)*(charL/beta)*sqrt(g)*(sin(sc)/sqrt(cos(sc)))*((tan(sc2)-tan(sc))/(tan(sc)-tan(sc1)));

    rint=2.5e-3;

    aint=1*rint;
    aend=50*rint;

    Gam=0.05877;
    q0=3.5e-6;%Q/2pi
    uint=2*q0/rint^2;

    % 相关参数
%     gammastar=(tan(sc2)-tan(sc))/(tan(sc)-tan(sc1));
%     gamma=(beta/charL)*(g*cos(sc))^(1/2);
%     hs=(q0*gammastar/(gamma*t))^(2/5);
%     us=(gamma/gammastar)*hs^(3/2);
    Res=uint*(rint)^(1/2)/nu;
    Fr=uint/(g*cos(sc)*rint)^(1/2);
    H=1/((t)*cos(sc)); %1/X_2
    V=1/(t);%1/A

    % 便捷计算代号
    C=5/4; %Velocity profile shape coeff.
    F=(1/Fr)^2;
    M=3/2;
    Y=Gam;%Bed stress coeff.
    R=Res;
    
    dy=zeros(5,1);
    dy(1)=1i*w*y(1)-V*y(1)-V*y(2)-1i*B*H*y(3)-y(4);
    dy(2)=y(4);
    dy(3)=y(5);
    dy(4)=R*(-M*Y*F*y(1)-1i*w*y(2)+2*cos(sc)*H+(1/R)*(H^2)*(B^2)*y(2)+Y*F*y(2)+1i*B*H*y(3)+C*y(4)+F*(1i*w*y(1)-V*y(1)-V*y(2)-1i*B*H*y(3)-y(4)));
    dy(5)=R*(1i*F*B*H*y(1)-M*Y*F*y(1)-1i*w*y(3)+(1/R)*(H^2)*(B^2)*y(3)+Y*F*y(3)+y(5));
end
% function dy = aODEf(t,y)
%     dy=zeros(1,1);
%     dy(1)=y(1);
% end
