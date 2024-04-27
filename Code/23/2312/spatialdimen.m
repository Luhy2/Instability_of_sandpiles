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
    R=1/R1;
    
    ThreeDre=[];% store k2 series in 3Dplot

    %N1=100;
    N2=400;
    N3=100;
    %AA=linspace(0,1,N1);
    BB=linspace(0,20,N2);
    CC=linspace(0,2,N3);
    
    for k=1:N2
        B=BB(k);
        x1=[];
        x2=[];
        x3=[];
        x4=[];
        x5=[];
        %for j=1:N1
        for j=1:N3
            %A=AA(j);
            C=CC(j);
            %Y=-Ys.*(1+(1/3).*(hs.^2).*((A)^2+(H^2)*(B)^2-0*H*A*B));
            %R=(1/R1)*((A)^2+(H^2)*(B)^2);
                %aaa=1;            
                %bbb=(2*A*1i+2*R+A*T*1i+2*F*Y).*1i;                
                %ccc=-(-(A^2)*F*1i+(A^2)*T*2i+(A^2)*1i-(R^2)*1i-(F^2)*(Y^2)*1i+3*A*R+3*A*F*Y+A*R*T-F*R*Y*2i-(B^2)*F*(H^2)*1i+A*F*(M+T)*Y+B*F*H*M*Y).*1i;
                %ddd=A*B*F*H*M*(T-1)*Y*1i+B*F*H*M*R*Y+(A^2)*F*(T+M)*Y*1i+A*F*M*R*Y-(B^2)*(F^2)*(H^2)*Y*1i+B*(F^2)*H*M*(Y^2)+2*A*F*R*Y+A*(B^2)*F*(H^2)*T-(B^2)*F*(H^2)*R*1i-A*(B^2)*F*(H^2)-(A^2)*(F^2)*Y*1i+(A^2)*R*T*1i+(A^2)*F*(Y-R)*1i+A*(F^2)*(Y^2)*(M+1)-(A^3)*T+(A^2)*R*1i+A*(R^2)+(A^3)*F;
                aaa=R^2;
                bbb=R*1i + R*T*1i - C*R^2 - F*R*1i;
                ccc=F - T - C*R*3i - C*R*T*1i + 2*F*R*Y + 2*B^2*H^2*R^2 + F*M*R*Y;
                ddd=C + C^2*R*2i - F^2*Y*1i - C*F + 2*C*T + F*Y*1i + F*M*Y*1i + F*T*Y*1i + B^2*H^2*R*1i - B^2*F*H^2*R*2i + B^2*H^2*R*T*1i - 2*C*F*R*Y - 2*B^2*C*H^2*R^2 + B*F*H*M*R*Y;
                eee=F^2*Y^2 - 2*C^2 - C^2*T + F^2*M*Y^2 - C*F*Y*3i + B^4*H^4*R^2 - B^2*F*H^2 - B^2*C*H^2*R*3i + B^2*F*H^2*T - C*F*M*Y*1i - C*F*T*Y*1i - B^2*C*H^2*R*T*1i + 2*B^2*F*H^2*R*Y - B*F*H*M*Y*1i + B*F*H*M*T*Y*1i + B^2*F*H^2*M*R*Y;
                fff=- B^4*C*H^4*R^2 - B^4*F*H^4*R*1i + M*B^3*F*H^3*R*Y + B^2*C^2*H^2*R*2i - 2*B^2*C*F*H^2*R*Y - B^2*C*F*H^2 - B^2*F^2*H^2*Y*1i - M*B*C*F*H*Y*1i + M*B*F^2*H*Y^2 + C^3 + C^2*F*Y*2i - C*F^2*Y^2;

            coeffs=[aaa,bbb,ccc,ddd,eee,fff];
            result=roots(coeffs);
            x1(j)=imag(result(1,1));
            x2(j)=imag(result(2,1));
            x3(j)=imag(result(3,1));
            x4(j)=imag(result(4,1));
            x5(j)=imag(result(5,1));
            rx3(j)=real(result(3,1));
            rx1(j)=real(result(1,1));
            rx2(j)=real(result(2,1));
            rx4(j)=real(result(4,1));
            rx5(j)=real(result(5,1));
        end
        OneDre(k,:)=-x5(2:end);
        ThreeDre(k,:)=x3;
%     figure()
%     subplot(5,2,1)
%     plot(CC,x1)
%     hold on
%     subplot(5,2,2)
%     plot(CC,rx1)
%     hold on
%     subplot(5,2,3)
%     plot(CC,x2)
%     hold on
%     subplot(5,2,4)
%     plot(CC,rx2)
%     hold on
%     subplot(5,2,5)
%     plot(CC,x3)
%     hold on
%     subplot(5,2,6)
%     plot(CC,rx3)
%     hold on
%     subplot(5,2,7)
%     plot(CC,x4)
%     hold on
%     subplot(5,2,8)
%     plot(CC,rx4)
%     hold on
%     subplot(5,2,9)
%     plot(CC,x5)
%     hold on
%     subplot(5,2,10)
%     plot(CC,rx5)
%     hold on
    
    %figure()
    
    %plot(CC,-x5)
    %hold on
    
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
        
    figure()
    %subplot(2,1,1)
    %surf(CC(2:end),BB.*H,OneDre)
    
    %subplot(2,1,2)
    [Graph,c]=contourf(CC(2:end),BB.*H,OneDre,[0,max(OneDre)],'EdgeColor','none');
    shading interp
    hold on
    c.LineWidth = 1;
    [Graphb,c]=contour(CC(2:end),BB.*H,OneDre,'LevelList',[0 1e-2],'LineColor','black','LineWidth',1,'ShowText','on');
    [Graphbb,cbb]=contour(CC(2:end),BB.*H,OneDre,'LevelList',[0.013 0.014],'LineColor','black','LineWidth',1,'ShowText','on');
    legendstr{i}=['Res=',num2str(Res),'   ','Fr=',num2str(Fr)];
    legend(legendstr)
    xlabel('\Omega')
    ylabel('K_2')
    title('H=0.0294')
    
     %colormap hot
%     subplot(3,1,3)
%     gl=imagesc(OneDre);
%     set(gca,'ydir','normal')
%     
%     hold on
    
end


