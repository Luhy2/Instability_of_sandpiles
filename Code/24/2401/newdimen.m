clc
clear
NUM=1;
%NU=linspace(40e-6,40e-6,NUM);
for i=1:NUM
    %nu=NU(i);
    Fr=1.02;
    Res=8.45;
    Gam=0.05877;
    
    T=5/4;
    F=(1/Fr)^2;
    M=3/2;
    H=1/34;
    Y=Gam;
    R1=8.45; %Re*
    
    ThreeDre=[];% store k2 series in 3Dplot

    N1=100;
    N2=100;
    AA=linspace(0,1,N1);
    BB=linspace(0,40,N2);
%     N1=20;
%     N2=1;
%     AA=linspace(0,1,N1);
%     BB=linspace(0,0,N2);
    %A=100;
    %B=1;
    %figure()
    for k=1:N2
        B=BB(k);
        x1=[];
        x2=[];
        x3=[];
        for j=1:N1
            A=AA(j);

            %Y=-Ys.*(1+(1/3).*(hs.^2).*((A)^2+(H^2)*(B)^2-0*H*A*B));
            R=(1/R1)*((A)^2+(H^2)*(B)^2);
                aaa=1;            
                bbb=(2*A*1i+2*R+A*T*1i+2*F*Y).*1i;                
                ccc=-(-(A^2)*F*1i+(A^2)*T*2i+(A^2)*1i-(R^2)*1i-(F^2)*(Y^2)*1i+3*A*R+3*A*F*Y+A*R*T-F*R*Y*2i-(B^2)*F*(H^2)*1i+A*F*(M+T)*Y+B*F*H*M*Y).*1i;
                ddd=A*B*F*H*M*(T-1)*Y*1i+B*F*H*M*R*Y+(A^2)*F*(T+M)*Y*1i+A*F*M*R*Y-(B^2)*(F^2)*(H^2)*Y*1i+B*(F^2)*H*M*(Y^2)+2*A*F*R*Y+A*(B^2)*F*(H^2)*T-(B^2)*F*(H^2)*R*1i-A*(B^2)*F*(H^2)-(A^2)*(F^2)*Y*1i+(A^2)*R*T*1i+(A^2)*F*(Y-R)*1i+A*(F^2)*(Y^2)*(M+1)-(A^3)*T+(A^2)*R*1i+A*(R^2)+(A^3)*F;


            coeffs=[aaa,bbb,ccc,ddd];
            result=roots(coeffs);
            x1(j)=imag(result(1,1));
            x2(j)=imag(result(2,1));
            x3(j)=imag(result(3,1));
            rx3(j)=real(result(3,1));
            rx1(j)=real(result(1,1));
        end
        OneDre(k,:)=x1;
        ThreeDre(k,:)=x3;
%     subplot(1,1,1)
%     plot(AA,x1)
%     hold on
%     subplot(3,1,2)
%     plot(AA,x2)
%     hold on
%     subplot(3,1,3)
%     plot(AA,x3)
%     hold on
    
    %subplot(2,1,2)
    %plot(AA,rx1./AA)
    %hold on
    %subplot(3,1,2)
    %plot(AA,x2)
    %hold on
    %subplot(3,1,3)
    %plot(AA,rx3./AA)
    %hold on
    
    end
    
    
    
    
%    figure()
    %subplot(2,2,1)
    %surf(AA,BB.*H,ThreeDre)
   
    
    %subplot(2,2,2)
    
    %[Graph,c]=contour(AA,BB.*H,ThreeDre,[0,max(ThreeDre)]);
    %hold on
    %c.LineWidth = 1;
    
    %legendstr{i}=['Res=',num2str(Res),'   ','Fr=',num2str(Fr)];
    %legend(legendstr)
    
    %colormap winter
    
    subplot(2,2,3)
    surf(AA,BB.*H,OneDre)
   
    
    subplot(2,2,4)
    
    [Graph,c]=contour(AA,BB.*H,OneDre,[0,max(OneDre)]);
    hold on
    c.LineWidth = 1;
    
    legendstr{i}=['Res=',num2str(Res),'   ','Fr=',num2str(Fr)];
    legend(legendstr)
    
    colormap winter

    
end


