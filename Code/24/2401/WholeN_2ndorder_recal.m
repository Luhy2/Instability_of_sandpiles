clc
clear
tic;
    %rint=2.5e-3;
    %default 1-41
    Aint=10;
    Aend=300;
    %hmax=0.4040;%限定网格点间距最大值
    %Num=ceil(1+(Aend-Aint)/hmax);
    Num=100;
    %Num=100;%离散点个数
    
    %h=(Aend-Aint)/(Num-1);
    
    A=linspace(Aint,Aend,Num);
    
    %a=A.*rint;
    h=(Aend-Aint)/(Num-1);
    
    %NBB=30;
    NBB=1;
    BB=linspace(1e-4,2.6,NBB);%k2域
    %BBB=BB.^2;

for op=1:NBB
    B=BB(op);%循环写k2，暂时先放一个数
    
    sc1=deg2rad(20.9);
    sc2=deg2rad(32.76);
    sc=deg2rad(29); %倾角

    %charL=0.825e-3;
    %beta=0.136;
    Fr=1.02;
    Res=8.45;
    
    charL=0.825e-3;
    beta=0.136;
    g=9.8;

    nu=(2/9)*(charL/beta)*sqrt(g)*(sin(sc)/sqrt(cos(sc)))*((tan(sc2)-tan(sc))/(tan(sc)-tan(sc1)));
    
    %Gam=0.05877; %Big Gamma
    gammastar=(tan(sc2)-tan(sc))/(tan(sc)-tan(sc1));
    gamma=(beta/charL)*(g*cos(sc))^(1/2);
    
    % Downstream info.
    usref=(Fr*Res*nu*(g*cos(sc))^(1/2))^(1/2);
    hsref=(Res/Fr)*nu/(g*cos(sc))^(1/2);
    Rint=hsref;
    aend=Aend*Rint;
    a=A.*Rint; %real a
    
    q0=aend*cos(sc)*hsref*usref;
    hs=(q0*gammastar./(gamma.*a)).^(2/5);
    us=(gamma/gammastar).*hs.^(3/2);%gammastar是文章里的gamma
    
    rint=hsref;
    uint=usref;
    
    Us=us./uint;
    Hs=hs./rint;
    
    %eta=(tan(sc2)-tan(sc1))/gammastar;
    
    G=gammastar;
    
    eta=tan(sc1)+(tan(sc2)-tan(sc1))/(1+G);
    Gamp=(G./((1+G).^2)).*(tan(sc2)-tan(sc1));
    
    % 便捷计算代号
    C=5/4; %Velocity profile shape coeff.
    F=(1/Fr)^2;
    Y=Gamp;%Bed stress coeff.
    R=Res;
    X2=1/cos(sc);
    %L=[];
    %各项系数设置成矩阵调用
    L(1,:)=-Us./(2.*h);
    L(2,:)=Us./A;
    L(3,:)=Hs./A;
    L(4,:)=1i.*B.*X2.*Hs./A;
    L(5,:)=Hs;
    L(6,:)=Us./(2.*h);
    L(7,:)=-F./(2.*h);
    L(8,:)=(1./R).*(Hs.^(1/2))./(2.*h);
    L(9,:)=-(3/2)*Y*F./Hs;
    
    L101=(2.*Us./A)+(Y*F./Us);
    L102=((1/R).*(Hs.^(1/2)).*((B.*X2).^2)./(A.^2));
    
    L(10,:)=(2.*Us./A)+(Y*F./Us)+((1/R).*(Hs.^(1/2)).*((B.*X2).^2)./(A.^2));
    L(11,:)=1i*B*X2.*Us./A;
    L(12,:)=C.*Us;
    L(13,:)=F./(2.*h);
    L(14,:)=-(1/R).*(Hs.^(1/2))./(2.*h);
    L(15,:)=(1/R).*(Hs.^(1/2))./(2.*h);
    L(16,:)=1i*F*B*X2./A;
    
    L171=(eta*F./Us);
    L172=((1/R).*(Hs.^(1/2)).*((B.*X2).^2)./(A.^2));
    
    L(17,:)=(eta*F./Us)+((1/R).*(Hs.^(1/2)).*((B.*X2).^2)./(A.^2));
    L(18,:)=Us;
    L(19,:)=-(1/R).*(Hs.^(1/2))./(2.*h);
    L(20,:)=-1./(2.*h);
    L(21,:)=-1;
    L(22,:)=1./(2.*h);
    L(23,:)=-1./(2.*h);
    L(24,:)=-1;
    L(25,:)=1./(2.*h);

%     L=zeros(25,Num);%only for test No.
%     for xp=1:25
%         L(xp,:)=xp;
%     end
    
    %初始化特征矩阵
    MATSIZE=(Num-2)*5;
    MAT=zeros(MATSIZE);
    RHS=zeros(MATSIZE);%右端特征矩阵
    
    KK(1,1)=2./h; %边界差分系数
    KK(1,2)=-1./(2.*h);
    KK(1,3)=4/3;
    KK(1,4)=-1/3;

%     KK(1,1)=1000;%only for test No.
%     KK(1,2)=2000;
%     KK(1,3)=3000;
%     KK(1,4)=4000;


    k=1;%记录每5行的切换
    g=0;%记录每5行的切换次数66698

    for i=1:MATSIZE
        if k==6
            k=1;
            g=g+1;
        end
        for j=1:MATSIZE
            if k==1
                if g==0 %左边值条件i=1代入i=2
                    MAT(i,1+5*g)=L(2,g+2);
                    MAT(i,2+5*g)=L(3,g+2);
                    MAT(i,3+5*g)=L(4,g+2);
                    MAT(i,4+5*g)=L(5,g+2);
                    
                    MAT(i,1+5*(g+1))=L(6,g+2);
                elseif g==Num-3 %右边值条件i=N代入i=N-1 -> 3表示前后各减1并更换5行的次数自带减1
                    MAT(i,1+5*(g-1))=L(1,g+2)+L(6,g+2)*KK(1,4);
                    
                    MAT(i,1+5*g)=L(2,g+2)+L(6,g+2)*KK(1,3);
                    MAT(i,2+5*g)=L(3,g+2);
                    MAT(i,3+5*g)=L(4,g+2);
                    MAT(i,4+5*g)=L(5,g+2);
                else
                    MAT(i,1+5*(g-1))=L(1,g+2);
                    
                    MAT(i,1+5*g)=L(2,g+2);
                    MAT(i,2+5*g)=L(3,g+2);
                    MAT(i,3+5*g)=L(4,g+2);
                    MAT(i,4+5*g)=L(5,g+2);
                    
                    MAT(i,1+5*(g+1))=L(6,g+2);
                end
            end
            if k==2
                if g==0
                    MAT(i,1+5*g)=L(9,g+2);
                    MAT(i,2+5*g)=L(10,g+2)+L(8,g+2)*KK(1,1);
                    MAT(i,3+5*g)=L(11,g+2);
                    MAT(i,4+5*g)=L(12,g+2);
                    
                    MAT(i,1+5*(g+1))=L(13,g+2);
                    MAT(i,2+5*(g+1))=L(8,g+2)*KK(1,2);
                    MAT(i,4+5*(g+1))=L(14,g+2);
                elseif g==Num-3
                    MAT(i,1+5*(g-1))=L(7,g+2)+L(13,g+2)*KK(1,4);
                    MAT(i,4+5*(g-1))=L(8,g+2);
                    
                    MAT(i,1+5*g)=L(9,g+2)+L(13,g+2)*KK(1,3);
                    MAT(i,2+5*g)=L(10,g+2);
                    MAT(i,3+5*g)=L(11,g+2);
                    MAT(i,4+5*g)=L(12,g+2);
                else
                    MAT(i,1+5*(g-1))=L(7,g+2);
                    MAT(i,4+5*(g-1))=L(8,g+2);
                    
                    MAT(i,1+5*g)=L(9,g+2);
                    MAT(i,2+5*g)=L(10,g+2);
                    MAT(i,3+5*g)=L(11,g+2);
                    MAT(i,4+5*g)=L(12,g+2);
                    
                    MAT(i,1+5*(g+1))=L(13,g+2);
                    MAT(i,4+5*(g+1))=L(14,g+2);
                end
            end
            if k==3
                if g==0
                    MAT(i,1+5*g)=L(16,g+2);
                    MAT(i,3+5*g)=L(17,g+2)+L(15,g+2)*KK(1,1);
                    MAT(i,5+5*g)=L(18,g+2);
                    
                    MAT(i,3+5*(g+1))=L(15,g+2)*KK(1,2);
                    MAT(i,5+5*(g+1))=L(19,g+2);                    
                elseif g==Num-3
                    MAT(i,5+5*(g-1))=L(15,g+2);
                    
                    MAT(i,1+5*g)=L(16,g+2);
                    MAT(i,3+5*g)=L(17,g+2);
                    MAT(i,5+5*g)=L(18,g+2);
                else
                    MAT(i,5+5*(g-1))=L(15,g+2);
                    
                    MAT(i,1+5*g)=L(16,g+2);
                    MAT(i,3+5*g)=L(17,g+2);
                    MAT(i,5+5*g)=L(18,g+2);
                    
                    MAT(i,5+5*(g+1))=L(19,g+2);
                end
            end
            if k==4
                if g==0
                    MAT(i,4+5*g)=L(21,g+2);
                    
                    MAT(i,2+5*(g+1))=L(22,g+2);
                elseif g==Num-3
                    MAT(i,2+5*(g-1))=L(20,g+2)+L(22,g+2)*KK(1,4);
                    
                    MAT(i,2+5*g)=L(22,g+2)*KK(1,3);
                    MAT(i,4+5*g)=L(21,g+2);
                else
                    MAT(i,2+5*(g-1))=L(20,g+2);
                    
                    MAT(i,4+5*g)=L(21,g+2);
                    
                    MAT(i,2+5*(g+1))=L(22,g+2);
                end
            end
            if k==5
                if g==0
                    MAT(i,5+5*g)=L(24,g+2);
                    
                    MAT(i,3+5*(g+1))=L(25,g+2);
                elseif g==Num-3
                    MAT(i,3+5*(g-1))=L(23,g+2)+L(25,g+2)*KK(1,4);
                    
                    MAT(i,3+5*g)=L(25,g+2)*KK(1,3);
                    MAT(i,5+5*g)=L(24,g+2);
                else
                    MAT(i,3+5*(g-1))=L(23,g+2);
                    
                    MAT(i,5+5*g)=L(24,g+2);
                    
                    MAT(i,3+5*(g+1))=L(25,g+2);
                end
            end
        end
        k=k+1;
    end
    opg=1;%同样是换5行计次
    for i=1:MATSIZE
        for j=1:MATSIZE
            if opg==6
                opg=1;
            end
            if i==j
                if opg<4
                    RHS(i,j)=1;
                else
                    RHS(i,j)=0;
                end
            else
                RHS(i,j)=0;
            end
        end
        opg=opg+1;
    end
    RRHS=1i.*RHS;
    QD=diag(RRHS);%用于检验右端特征矩阵的对角元素
    
    [V,D]=eig(MAT,RRHS);
    [VNEW,DNEW]=eigs(MAT,RRHS,5,'smallestabs');
    WWNEW=diag(DNEW);
    REE=diag(D);
    RE=REE.';
    matsize=MATSIZE-2*(Num-2);%减去不合理特征根由降阶产生的每个内点2个：排列在REE最后
    for pp=1:matsize
        WW(pp,op)=RE(1,pp);
        
    end
    Vnew=V(:,1:matsize);
    
%   
%     REEE=[];
%     uppk=1;
%     for i=1:MATSIZE-2*(Num-2)
%         
%         if real(REE(i,1))>1e-4
%             if real(REE(i,1))<1
%                 REEE(uppk,1)=REE(i,1);
%                 uppk=uppk+1;
%             end
%         end
%     end
%     for pp=1:uppk-1
%         WW(pp,op)=REEE(pp,1);
%     end
%
    %REE=[];
    
end

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
%NUUM=num2str(Num);

% for j=m:m
%     subplot(1,1,1)
%     plot(BB.*rint,imag(SWW(j,:)),'-o')
%     hold on
%     xlabel('K_2')
%     ylabel('\omega_I')
%     %legend('N=41','N=61','N=81')
% %     plot(BBB.*rint,imag(SWW(j,:)),'-o')
% %     hold on
% end

elapsedTime = toc;
disp(['程序运行时间：' num2str(elapsedTime) '秒']);

%%
%figure()
realWW=real(WW);
imagWW=imag(WW);
plot(realWW,imagWW,'o')
hold on

% plot(real(WWNEW),imag(WWNEW),'o')
% hold on
% for i = 1:1:matsize
%     text(realWW(i), imagWW(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
% end
% 
% hold off;
%%
clc
    
    figure()
    %clc
    TRY=200;
    for i=1:MATSIZE/5
        AlphaH(i,1)=Vnew((i-1)*5+1,TRY);
        AlphaU(i,1)=Vnew((i-1)*5+2,TRY);
        AlphaV(i,1)=Vnew((i-1)*5+3,TRY);
        PU(i,1)=Vnew((i-1)*5+4,TRY);
        PV(i,1)=Vnew((i-1)*5+5,TRY);
    end
    
    subplot(5,1,1)
    plot(A(2:end-1),real(AlphaH),'-o')
    hold on
    plot(A(2:end-1),imag(AlphaH))
    hold on

    subplot(5,1,2)
    plot(A(2:end-1),real(AlphaU),'-o')
    hold on
    plot(A(2:end-1),imag(AlphaU))
    hold on

    subplot(5,1,3)
    plot(A(2:end-1),real(AlphaV),'-o')
    hold on
    plot(A(2:end-1),imag(AlphaV))
    hold on


    subplot(5,1,4)
    plot(A(2:end-1),real(PU),'-o')
    hold on
    plot(A(2:end-1),imag(PU))
    hold on

    subplot(5,1,5)
    plot(A(2:end-1),real(PV),'-o')
    hold on
    plot(A(2:end-1),imag(PV))
    hold on
    
    
    
    
    

