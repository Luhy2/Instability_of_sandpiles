clc
clear
% Origin input parameters
Fr=1.02;
Res=8.45;
Gam=0.05877;
A=300; %nondimension length by rint
sc1=deg2rad(20.9);
sc2=deg2rad(32.76);
sc=deg2rad(29); %倾角
charL=0.825e-3;
beta=0.136;
g=9.8;
C=5/4; %profile shape coeff.
    nu=(2/9)*(charL/beta)*sqrt(g)*(sin(sc)/sqrt(cos(sc)))*((tan(sc2)-tan(sc))/(tan(sc)-tan(sc1)));
    gammastar=(tan(sc2)-tan(sc))/(tan(sc)-tan(sc1));
    gamma=(beta/charL)*(g*cos(sc))^(1/2);

% Downstream info.
    us=(Fr*Res*nu*(g*cos(sc))^(1/2))^(1/2);
    hs=(Res/Fr)*nu/(g*cos(sc))^(1/2);
    
% Incoming flow char.
Rint=2.5e-2;
a=A*Rint; %real a
%q0=3.5e-7; %Q/2pi here determined by init Fr and Res from 'infty' downstream
q0=a*cos(sc)*hs*us;

%uint=2*q0/(rint^2); %only ref velocity
rint=hs;
uint=us;


%New nondimen Numbers
F0=uint/(g*rint*cos(sc))^(1/2);
R0=(uint*(rint)^(1/2))/nu;

Us=us/uint;
Hs=hs/rint;

    %eta=(tan(sc2)-tan(sc1))/gammastar;
    eta=tan(sc1)+(tan(sc2)-tan(sc1))/(1+gammastar);
    G=gammastar;
    Gamp=(G./((1+G).^2)).*(tan(sc2)-tan(sc1));
%     eta = Gamp;
% Eigenvalue sol.s
N1=100;
N2=100;
AA=linspace(0,0.5,N1);
% BB=linspace(-5,10,N2);
BB = linspace(0,0,N2);
for k=1:N2
    B=BB(k);
    x1=[];
    x2=[];
    x3=[];
    for j=1:N1
        Ak=AA(j);
        % Only for test

        L(1,1)=1i*Ak*Us+Us/A;
        L(1,2)=1i*Ak*Hs+Hs/A;
%         L(1,1)=1i*Ak*Us;
%         L(1,2)=1i*Ak*Hs;
        L(1,3)=1i*B*Hs/(A*cos(sc));
        L(2,1)=1i*Ak*(1/F0^2)-(3/2)*(Gamp/F0^2)/Hs;
        L(2,2)=1i*C*Ak*Us+2*Us/A+(Hs^(1/2))*(1/R0)*(Ak^2+(B/(A*cos(sc)))^2)+(Gamp/F0^2)/Us;
%         L(2,2)=1i*C*Ak*Us+(Hs^(1/2))*(1/R0)*(Ak^2+(B/(A*cos(sc)))^2)+(Gamp/F0^2)/Us;
        L(2,3)=1i*B*Us/(A*cos(sc));

        L(3,1)=1i*(B/F0^2)/(A*cos(sc));
%         L(3,1)=1i*(B/F0^2)/(A*cos(sc))-(3/2)*(Gamp/F0^2)/Hs;
        L(3,2)=0;
        L(3,3)=1i*Ak*Us+(Hs^(1/2))*(1/R0)*(Ak^2+(B/(A*cos(sc)))^2)+(eta/F0^2)/Us;
%         L(3,3)=1i*Ak*Us+(Hs^(1/2))*(1/R0)*(Ak^2+(B/(A*cos(sc)))^2)+(Gamp/F0^2)/Us; 

        RHS=1i.*eye(3);
        
        [V,D]=eig(L,RHS);
        RD=diag(D);
        
        WW(:,j)=RD;
        
    end
%end

% 获取矩阵的大小
[m, n] = size(WW);
% 
% 初始化排序后的复矩阵
SWW = zeros(m, n);
% 逐列排序
for col = 1:n
    % 提取当前列的虚部
    imag_parts = imag(WW(:, col));
    % 获取排序后的索引
    [~, sorted_indices] = sort(imag_parts);
    % 根据排序后的索引重新排列当前列
    SWW(:, col) = WW(sorted_indices, col);
end

%Draw
    for j=m:m
        subplot(2,1,1)
        plot(AA,imag(SWW(j,:)),'-o')
        hold on
        xlabel('K_1')
        ylabel('\omega_I')
        subplot(2,1,2)
        plot(AA(2:end),real(SWW(j,(2:end))./AA(2:end)),'-o')
        hold on
        %legend('N=41','N=61','N=81')
    %     plot(BBB.*rint,imag(SWW(j,:)),'-o')
    %     hold on
    end

    OneDre(k,:)=imag(SWW(m,:));
end

%     subplot(2,1,1)
%     surf(AA,BB*(1/(A*cos(sc))),OneDre)
%    
%     
%     subplot(2,1,2)
%     [Graph,c]=contour(AA,BB*(1/(A*cos(sc))),OneDre,[0,max(OneDre)]);
%     %[Graph,c]=contourf(AA,BB*(1/(A*cos(sc))),OneDre,[0,max(OneDre)],'EdgeColor','none');
%     shading interp
%     hold on
%     c.LineWidth = 1;
%     [Graphb,c]=contour(AA,BB*(1/(A*cos(sc))),OneDre,'LevelList',[0.0088 0.0092],'LineColor','black','LineWidth',1,'ShowText','on');
%     [Graphbb,cbb]=contour(AA,BB*(1/(A*cos(sc))),OneDre,'LevelList',[0.0080,0.0088],'LineColor','black','LineWidth',1,'ShowText','on');
%     [Graphbbb,cbbb]=contour(AA,BB*(1/(A*cos(sc))),OneDre,'LevelList',[0,0],'LineColor','black','LineWidth',1,'ShowText','on');
%     legendstr=['Res=',num2str(Res),'   ','Fr=',num2str(Fr)];
%     legend(legendstr)
%     xlabel('K_1')
%     ylabel('K_2')
%     colormap winter

    
    
%     [Graph,c]=contourf(CC(2:end),BB.*H,OneDre,[0,max(OneDre)],'EdgeColor','none');
%     shading interp
%     hold on
%     c.LineWidth = 1;
%     [Graphb,c]=contour(CC(2:end),BB.*H,OneDre,'LevelList',[0 1e-2],'LineColor','black','LineWidth',1,'ShowText','on');
%     [Graphbb,cbb]=contour(CC(2:end),BB.*H,OneDre,'LevelList',[0.013 0.014],'LineColor','black','LineWidth',1,'ShowText','on');
%     legendstr{i}=['Res=',num2str(Res),'   ','Fr=',num2str(Fr)];
%     legend(legendstr)
%     xlabel('\Omega')
%     ylabel('K_2')
%     title('1/X_2=',H)

