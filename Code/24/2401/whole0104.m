clc
clear
    rint=2.5e-3;
    Aint=1;
    Aend=41;
    Num=5;
    A=linspace(Aint,Aend,Num);
    a=A.*rint;
    h=(Aend-Aint)/(Num-1);
    
    NBB=100;
    BB=linspace(0,5,NBB);%k2域
for op=1:NBB
    B=BB(op);%循环写k2，暂时先放一个数
    
    sc1=deg2rad(20.9);
    sc2=deg2rad(32.76);
    sc=deg2rad(29); %倾角
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
    
    %各项系数
    L1=(Us./h)+(Us./A);
    L2=-Us./h;
    L3=(Hs./h)+(Hs./A);
    L4=-Hs./h;
    L5=1i.*B.*X2.*Hs./A;
    L6=(F/h)-(3/2).*Y.*F./Hs;
    L7=-F/h;
    L8=-(1/(R*h^2)).*(Hs).^(1/2);
    L9=((C/h).*Us+(2./A).*Us+(2/(R*h^2)).*(Hs.^(1/2))+(1/R).*((B.*X2).^2).*(1./A.^2).*Hs.^(1/2)+Gamp.*F./Us);
    %L9=((C/h).*Us+(2/(R*h^2)).*(Hs.^(1/2))+(1/R).*((B.*X2).^2).*(1./A.^2).*Hs.^(1/2)+Gamp.*F./Us);
    L10=-(C/h).*Us-(1/(R*h^2)).*Hs.^(1/2);
    L11=1i.*B.*X2.*Us./A;
    %L11=[0 0 0 0 0];
    L12=1i.*F.*B.*X2./A;
    L13=-(1/(R*h^2)).*Hs.^(1/2);
    L14=(Us./h)+((2/(R*h^2)).*Hs.^(1/2))+(1/R).*((B.*X2).^2).*(1./A.^2).*Hs.^(1/2)+eta.*F./Us;
    L15=-Us./h-(1/(R*h^2)).*Hs.^(1/2);
    
    %初始化特征矩阵
    MAT=zeros((Num-2)*3);
    
    %暂时手填
    MAT(1,1)=L1(1,2);
    MAT(1,2)=L3(1,2);
    MAT(1,3)=L5(1,2);
    
    MAT(2,1)=L6(1,2);
    MAT(2,2)=L9(1,2);
    MAT(2,3)=L11(1,2);
    MAT(2,5)=L8(1,2);
    
    MAT(3,1)=L12(1,2);
    MAT(3,3)=L14(1,2);
    MAT(3,6)=L13(1,2);
    
    MAT(4,1)=L2(1,3);
    MAT(4,2)=L4(1,3);
    MAT(4,4)=L1(1,3);
    MAT(4,5)=L3(1,3);
    MAT(4,6)=L5(1,3);
    
    MAT(5,1)=L7;
    MAT(5,2)=L10(1,3);
    MAT(5,4)=L6(1,3);
    MAT(5,5)=L9(1,3);
    MAT(5,6)=L11(1,3);
    MAT(5,8)=L8(1,3);
    
    MAT(6,3)=L15(1,3);
    MAT(6,4)=L12(1,3);
    MAT(6,6)=L14(1,3);
    MAT(6,9)=L13(1,3);
    
    MAT(7,4)=L2(1,4);
    MAT(7,5)=L4(1,4);
    MAT(7,7)=L1(1,4);
    MAT(7,8)=L3(1,4);
    MAT(7,9)=L5(1,4);
    
    MAT(8,4)=L7;
    MAT(8,5)=L10(1,4);
    MAT(8,7)=L6(1,4);
    MAT(8,8)=L8(1,4)+L9(1,4);
    MAT(8,9)=L11(1,4);
    
    MAT(9,6)=L15(1,4);
    MAT(9,7)=L12(1,4);
    MAT(9,9)=L13(1,4)+L14(1,4);
    
    
    RHS=1i.*eye((Num-2)*3);
    
    [V,D]=eig(MAT,RHS);
    
    REE=diag(D);
    RE=REE.';
    
    WW(1,op)=RE(1,1);
    WW(2,op)=RE(1,2);
    WW(3,op)=RE(1,3);
    WW(4,op)=RE(1,4);
    WW(5,op)=RE(1,5);
    WW(6,op)=RE(1,6);
    WW(7,op)=RE(1,7);
    WW(8,op)=RE(1,8);
    WW(9,op)=RE(1,9);
    
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
figure()
plot(BB.*rint,imag(SWW(9,:)),'o')
hold on
plot(BB.*rint,imag(SWW(8,:)),'o')
hold on
plot(BB.*rint,imag(SWW(7,:)),'o')
hold on
plot(BB.*rint,imag(SWW(6,:)),'o')
hold on
plot(BB.*rint,imag(SWW(5,:)),'o')
hold on
plot(BB.*rint,imag(SWW(4,:)),'o')
hold on
plot(BB.*rint,imag(SWW(3,:)),'o')
hold on
plot(BB.*rint,imag(SWW(2,:)),'o')
hold on
plot(BB.*rint,imag(SWW(1,:)),'o')
hold on
legend('1','2','3','4','5','6','7','8','9')


figure()
plot(BB.*rint,real(SWW(9,:)),'o')
hold on
plot(BB.*rint,real(SWW(8,:)),'o')
hold on
plot(BB.*rint,real(SWW(7,:)),'o')
hold on
plot(BB.*rint,real(SWW(6,:)),'o')
hold on
plot(BB.*rint,real(SWW(5,:)),'o')
hold on
plot(BB.*rint,real(SWW(4,:)),'o')
hold on
plot(BB.*rint,real(SWW(3,:)),'o')
hold on
plot(BB.*rint,real(SWW(2,:)),'o')
hold on
plot(BB.*rint,real(SWW(1,:)),'o')
hold on
legend('1','2','3','4','5','6','7','8','9')
    
% figure()
% plot(BB.*rint,imag(WW(9,:)),'o')
% hold on
% plot(BB.*rint,imag(WW(8,:)),'o')
% hold on
% plot(BB.*rint,imag(WW(7,:)),'o')
% hold on
% plot(BB.*rint,imag(WW(6,:)),'o')
% hold on
% plot(BB.*rint,imag(WW(5,:)),'o')
% hold on
% plot(BB.*rint,imag(WW(4,:)),'o')
% hold on
% plot(BB.*rint,imag(WW(3,:)),'o')
% hold on
% plot(BB.*rint,imag(WW(2,:)),'o')
% hold on
% plot(BB.*rint,imag(WW(1,:)),'o')
% hold on
% legend('1','2','3','4','5','6','7','8','9')
% 
% 
% figure()
% plot(BB.*rint,real(WW(9,:)),'o')
% hold on
% plot(BB.*rint,real(WW(8,:)),'o')
% hold on
% plot(BB.*rint,real(WW(7,:)),'o')
% hold on
% plot(BB.*rint,real(WW(6,:)),'o')
% hold on
% plot(BB.*rint,real(WW(5,:)),'o')
% hold on
% plot(BB.*rint,real(WW(4,:)),'o')
% hold on
% plot(BB.*rint,real(WW(3,:)),'o')
% hold on
% plot(BB.*rint,real(WW(2,:)),'o')
% hold on
% plot(BB.*rint,real(WW(1,:)),'o')
% hold on
% legend('1','2','3','4','5','6','7','8','9')
    
    
    
    
    
    
    