clc
clear
    rint=2.5e-3;
    Aint=1;
    Aend=41;
    Num=31;
    A=linspace(Aint,Aend,Num);
    a=A.*rint;
    h=(Aend-Aint)/(Num-1);
    
    NBB=100;
    BB=linspace(0,5,NBB);%k2域

for op=1:NBB
    B=BB(op);%循环写k2，暂时先放一个数
    
    sc1=deg2rad(20.9);
    sc2=deg2rad(32.76);
    sc=deg2rad(25); %倾角

    charL=0.825e-3;
    beta=0.136;
    g=9.8;

    nu=(2/9)*(charL/beta)*sqrt(g)*(sin(sc)/sqrt(cos(sc)))*((tan(sc2)-tan(sc))/(tan(sc)-tan(sc1)));

    %Gam=0.05877; %Big Gamma
    q0=3.5e-6;%Q/2pi
    uint=2*q0/rint^2;
    
    gammastar=(tan(sc2)-tan(sc))/(tan(sc)-tan(sc1));
    gamma=(beta/charL)*(g*cos(sc))^(1/2);
    %暂时用参数化表达hs与us
    hs=(q0*gammastar./(gamma.*a)).^(2/5);
    us=(gamma/gammastar).*hs.^(3/2);%gammastar是文章里的gamma
    Frlocal=us./(g.*hs.*cos(sc)).^(1/2);
    Hs=hs./rint;
    Us=us./uint;
    
    Res=uint*(rint)^(1/2)/nu;
    Fr=uint/(g*cos(sc)*rint)^(1/2);
    
    eta=(tan(sc2)-tan(sc1))/gammastar;
    
    %eta=0;
    %G=gammastar.*(Hs.^(3/2))./Us;
    G=gammastar;
    Gamp=(G./((1+G).^2)).*(tan(sc2)-tan(sc1));
    
    % 便捷计算代号
    C=5/4; %Velocity profile shape coeff.
    F=(1/Fr)^2;
    Y=Gamp;%Bed stress coeff.
    R=Res;
    X2=1/cos(sc);
    %L=[];
    %各项系数设置成矩阵调用
    L(1,:)=(Us./h)+(Us./A);
    L(2,:)=-Us./h;
    L(3,:)=(Hs./h)+(Hs./A);
    L(4,:)=-Hs./h;
    L(5,:)=1i.*B.*X2.*Hs./A;
    L(6,:)=(F/h)-(3/2).*Y.*F./Hs;
    L(7,:)=-F/h;
    L(8,:)=-(1/(R*h^2)).*(Hs).^(1/2);
    L(9,:)=((C/h).*Us+(2./A).*Us+(2/(R*h^2)).*(Hs.^(1/2))+(1/R).*((B.*X2).^2).*(1./A.^2).*Hs.^(1/2)+Gamp.*F./Us);
    %L(9,:)=((C/h).*Us+(2/(R*h^2)).*(Hs.^(1/2))+(1/R).*((B.*X2).^2).*(1./A.^2).*Hs.^(1/2)+Gamp.*F./Us);
    L(10,:)=-(C/h).*Us-(1/(R*h^2)).*Hs.^(1/2);
    L(11,:)=1i.*B.*X2.*Us./A;
    %L(11,:)=zeros(1,Num);
    L(12,:)=1i.*F.*B.*X2./A;
    L(13,:)=-(1/(R*h^2)).*Hs.^(1/2);
    L(14,:)=(Us./h)+((2/(R*h^2)).*Hs.^(1/2))+(1/R).*((B.*X2).^2).*(1./A.^2).*Hs.^(1/2)+eta.*F./Us;
    L(15,:)=-Us./h-(1/(R*h^2)).*Hs.^(1/2);
    
    %初始化特征矩阵
    MATSIZE=(Num-2)*3;
    MAT=zeros(MATSIZE);
    k=1;%记录每3行的切换
    g=0;%记录每3行的切换次数
    for i=1:MATSIZE
        if k==4
            k=1;
            g=g+1;
        end
        for j=1:MATSIZE
            if k==1
                if g==0

                    MAT(i,1+3*g)=L(1,g+2);
                    MAT(i,2+3*g)=L(3,g+2);
                    MAT(i,3+3*g)=L(5,g+2);
                else
                    MAT(i,1+3*(g-1))=L(2,g+2);
                    MAT(i,2+3*(g-1))=L(4,g+2);
                    MAT(i,1+3*g)=L(1,g+2);
                    MAT(i,2+3*g)=L(3,g+2);
                    MAT(i,3+3*g)=L(5,g+2);
                end
            end
            if k==2
                if g==0
                    
                    MAT(i,1+3*g)=L(6,g+2);
                    MAT(i,2+3*g)=L(9,g+2);
                    MAT(i,3+3*g)=L(11,g+2);
                    MAT(i,5+3*g)=L(8,g+2);
                elseif g==Num-3
                    MAT(i,1+3*(g-1))=L(7,g+2);
                    MAT(i,2+3*(g-1))=L(10,g+2);
                    MAT(i,1+3*g)=L(6,g+2);
                    MAT(i,2+3*g)=L(8,g+2)+L(9,g+2);
                    MAT(i,3+3*g)=L(11,g+2);
                else
                    MAT(i,1+3*(g-1))=L(7,g+2);
                    MAT(i,2+3*(g-1))=L(10,g+2);
                    MAT(i,1+3*g)=L(6,g+2);
                    MAT(i,2+3*g)=L(9,g+2);
                    MAT(i,3+3*g)=L(11,g+2);
                    MAT(i,5+3*g)=L(8,g+2);
                end
            end
            if k==3
                if g==0
                    
                    MAT(i,1+3*g)=L(12,g+2);
                    MAT(i,3+3*g)=L(14,g+2);
                    MAT(i,6+3*g)=L(13,g+2);
                elseif g==Num-3
                    MAT(i,3+3*(g-1))=L(15,g+2);
                    MAT(i,1+3*g)=L(12,g+2);
                    MAT(i,3+3*g)=L(13,g+2)+L(14,g+2);
                else
                    MAT(i,3+3*(g-1))=L(15,g+2);
                    MAT(i,1+3*g)=L(12,g+2);
                    MAT(i,3+3*g)=L(14,g+2);
                    MAT(i,6+3*g)=L(13,g+2);
                end
            end
            
        end
        k=k+1;
    end

    RHS=1i.*eye((Num-2)*3);
    
    [V,D]=eig(MAT,RHS);
    
    REE=diag(D);
    RE=REE.';
    for pp=1:MATSIZE
        WW(pp,op)=RE(1,pp);
    end
end

% 获取矩阵的大小
[m, n] = size(WW);
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

subplot(2,1,1)
plot(BB.*rint,imag(SWW(MATSIZE,:)),'o')
hold on

% plot(BB.*rint,imag(SWW(MATSIZE-20,:)),'o')
% hold on
% 
% plot(BB.*rint,imag(SWW(MATSIZE-50,:)),'o')
% hold on
xlabel('K_2')
ylabel('\omega_I')
legend('N=30','N=100','N=200')

subplot(2,1,2)
plot(BB.*rint,real(SWW(MATSIZE,:)),'o')
hold on

% plot(BB.*rint,real(SWW(MATSIZE-20,:)),'o')
% hold on
% 
% plot(BB.*rint,real(SWW(MATSIZE-50,:)),'o')
% hold on
xlabel('K_2')
ylabel('\omega_R')
legend('N=30','N=100','N=200')
% subplot(3,1,2)
% plot(BB,real(SWW(MATSIZE,:)),'o')
% hold on
% subplot(3,1,3)
% plot(BB,real(SWW(MATSIZE,:))./BB,'o')
% hold on
% 
% 
% figure()
% plot(BB.*rint,real(SWW(9,:)),'o')
% hold on
% plot(BB.*rint,real(SWW(8,:)),'o')
% hold on
% plot(BB.*rint,real(SWW(7,:)),'o')
% hold on
% plot(BB.*rint,real(SWW(6,:)),'o')
% hold on
% plot(BB.*rint,real(SWW(5,:)),'o')
% hold on
% plot(BB.*rint,real(SWW(4,:)),'o')
% hold on
% plot(BB.*rint,real(SWW(3,:)),'o')
% hold on
% plot(BB.*rint,real(SWW(2,:)),'o')
% hold on
% plot(BB.*rint,real(SWW(1,:)),'o')
% hold on
% legend('1','2','3','4','5','6','7','8','9')
%     
% % figure()
% % plot(BB.*rint,imag(WW(9,:)),'o')
% % hold on
% % plot(BB.*rint,imag(WW(8,:)),'o')
% % hold on
% % plot(BB.*rint,imag(WW(7,:)),'o')
% % hold on
% % plot(BB.*rint,imag(WW(6,:)),'o')
% % hold on
% % plot(BB.*rint,imag(WW(5,:)),'o')
% % hold on
% % plot(BB.*rint,imag(WW(4,:)),'o')
% % hold on
% % plot(BB.*rint,imag(WW(3,:)),'o')
% % hold on
% % plot(BB.*rint,imag(WW(2,:)),'o')
% % hold on
% % plot(BB.*rint,imag(WW(1,:)),'o')
% % hold on
% % legend('1','2','3','4','5','6','7','8','9')
% % 
% % 
% % figure()
% % plot(BB.*rint,real(WW(9,:)),'o')
% % hold on
% % plot(BB.*rint,real(WW(8,:)),'o')
% % hold on
% % plot(BB.*rint,real(WW(7,:)),'o')
% % hold on
% % plot(BB.*rint,real(WW(6,:)),'o')
% % hold on
% % plot(BB.*rint,real(WW(5,:)),'o')
% % hold on
% % plot(BB.*rint,real(WW(4,:)),'o')
% % hold on
% % plot(BB.*rint,real(WW(3,:)),'o')
% % hold on
% % plot(BB.*rint,real(WW(2,:)),'o')
% % hold on
% % plot(BB.*rint,real(WW(1,:)),'o')
% % hold on
% % legend('1','2','3','4','5','6','7','8','9')
    
 
    
    
    
    
    