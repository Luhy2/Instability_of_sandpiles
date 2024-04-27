clc
clear

N=1400;

xint=0;
xend=2;

x=linspace(xint,xend,N);
h=(xend-xint)/(N-1);

K1=1/h^2;
K2=-2/(h^2);
K3=1/h^2;

% 右边界一阶导为零
R1=3/(2*h);
R2=-2/h;
R3=1/(2*h);

% 右边界值为零
% R1=1;
% R2=0;
% R3=0;

LHS=zeros(N);
rHS=eye(N);
RHS=-1.*rHS;
RHS(1,1)=0;
RHS(N,N)=0;

for i=1:N
    if i==1
        LHS(i,1)=1;
    elseif i==N
        LHS(i,N)=R1;
        LHS(i,N-1)=R2;
        LHS(i,N-2)=R3;
    else
        LHS(i,i-1)=K3;
        LHS(i,i)=K2;
        LHS(i,i+1)=K1;
    end
end

% 边界值为零条件下的子矩阵
% SubLHS=LHS(2:end-1,2:end-1);
% SubRHS=RHS(2:end-1,2:end-1);

% 边界一阶导为零条件下的子矩阵

SubLHS=LHS(2:end-1,2:end-1);
SubLHS(end,end)=K2-K1.*(R2./R1);
SubLHS(end,end-1)=K3-K1.*(R3./R1);

SubRHS=RHS(2:end-1,2:end-1);


% [V,D]=eig(LHS,RHS);
[V,D]=eigs(LHS,RHS,N,'smallestabs');
[SubV,SubD]=eigs(SubLHS,SubRHS,N-2,'smallestabs');
Omegalist=diag(D);
SubOmegalist=diag(SubD);
% realOme=real(Omegalist);
% imagOme=imag(Omegalist);

Length=length(Omegalist);
SubLen=length(SubOmegalist);
LST=linspace(1,Length,Length);
SubLST=linspace(1,SubLen,SubLen);
%%
% verify omegalist and \sqrt(k)*x_N=(m+0.5)\pi

% figure()
plot(LST,Omegalist(1:end,1),'o')
hold on

plot(SubLST,SubOmegalist(1:end,1),'^')
hold on

sqrtk=((LST+0.5).*pi)./xend;
sqrtk=((LST).*pi)./xend;
doubk=sqrtk.^2;
plot(LST,doubk,'o')
hold on

%%
clc
figure()
TRY=469;
%  plot(x,real(V(:,TRY)),'-o')
%  hold on
%  plot(x,-13.23*V(:,TRY),'-o')
sqrtkTRY=((TRY).*pi)./xend;

plot(x,V(:,TRY),'-o')
 hold on 
plot(x(1,2:end-1),SubV(:,TRY),'^')
hold on


% xc=linspace(xint,xend,20*N);
% doubk=sqrtkTRY.^2;
% COMP1=sin((sqrtkTRY).*xc);
% COMP2=cos((sqrtkTRY).*xc);
% plot(xc,COMP1,'o')
% hold on
%%
 xc=linspace(xint,xend,20*N);
 
 COMP=1*sin(xc.*sqrt(Omegalist(TRY,1)));
 COMP2=-1*sin(xc.*sqrtkTRY);

 plot(xc,COMP2,'LineWidth', 1)
 hold on
 title('sqrtkTRY=',sqrtkTRY)
 
 
 legend(['Eig=',num2str(sqrt(Omegalist(TRY,1)))])




