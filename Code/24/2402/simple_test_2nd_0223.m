clc
clear

% N=400;

N=100;

% xint=0;
% xend=(5/2).*pi;
% xint=-1;
% xend=1;
xint=0;
xend=20;

x=linspace(xint,xend,N);
h=(xend-xint)/(N-1);

% K1=1/h^2;
% K2=1-(2/(h^2));
% K3=1/h^2;

K1=1/h^2;
K2=-2/(h^2);
K3=1/h^2;


R1=3/(2*h);
R2=-2/h;
R3=1/(2*h);

% R1=3;
% R2=-4;
% R3=1;

% R1=1/h;
% R2=-1/h;
% R3=0;

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
%LLHS=rref(LHS);

% 三对角预处理子
P1 = LHS-tril(LHS,-2)-triu(LHS,2); 
% 下三角预处理子
P2 = tril(LHS,0);


% [V,D]=eig(LHS,RHS);
% [V,D]=eigs(LHS,RHS,N,'smallestabs');

[V,D]=eigs(P2\LHS,P2\RHS,N,'smallestabs');

% [V,D]=eigs(LHS,RHS,N);
Omegalist=diag(D);      
% realOme=real(Omegalist);
% imagOme=imag(Omegalist);
Length=length(Omegalist);
LST=linspace(1,Length,Length);

%%
% reduction method
LHSr=zeros(2*N);
rHSr=zeros(2*N);

M1=3/(2*h);
M2=-1;
M3=4/(2*h);
M4=1/(2*h);

Rr1=1/(2*h);
Rr2=-1;
Rr3=-1/(2*h);

for i=1:2*N
    if i==1
        LHSr(i,1)=1;
    elseif i==2
        LHSr(i,1)=-M1;
        LHSr(i,2)=M2;
        LHSr(i,3)=M3;
        LHSr(i,5)=-M4;
    elseif i==2*N
%         LHSr(i,i)=1;
    elseif i==2*N-1
        LHSr(i,2*N)=M1;
        LHSr(i,2*N-2)=-M3;
        LHSr(i,2*N-4)=M4;
        rHSr(i,i)=-1;
    else
        if mod(i,2)==1
            %odd
            LHSr(i,i-1)=Rr3;
            LHSr(i,i+3)=Rr1;
            rHSr(i,i)=-1;
        end
        if mod(i,2)==0
            %even
            LHSr(i,i)=Rr2;
            LHSr(i,i+1)=Rr1;
            LHSr(i,i-3)=Rr3;
        end
    end
end
LHSr(2*N,2*N-1)=1;
[Vr,Dr]=eigs(LHSr,rHSr,2*N,'smallestabs');
OmegalistR=diag(Dr);

%%
% reduction method diagrams
clc
LengthR=length(OmegalistR);
LSTr=linspace(1,LengthR,LengthR);

plot(LSTr,OmegalistR(1:end,1),'o')
hold on
%%
%figure()
        
TRY=6;
Vr1 = [];
Vr2 = [];


odd_indices = 1:2:length(Vr(:,TRY));
Vr1 = Vr(odd_indices,TRY);

even_indices = 2:2:length(Vr(:,TRY));
Vr2 = Vr(even_indices,TRY);

plot(x,Vr1,'o')
 hold on
%  plot(x,VR2,'o')
%  hold on 


%%
% verify omegalist and \sqrt(k)*x_N=(m+0.5)\pi


plot(LST,Omegalist(1:end,1),'o')
hold on
% plot(LST,OmegalistR(1:Length,1),'^')
% hold on

% sqrtk=((LST+0.5).*pi)./xend;
% sqrtk=((LST+0.5).*pi)./2;
% doubk=sqrtk.^2;
% plot(LST,doubk,'o')
% hold on
% legend('1000','400','200','Theoretical')
%%
clc
%figure()
TRY=3;
%  plot(x,real(V(:,TRY)),'-o')
%  hold on
%  plot(x,-13.23*V(:,TRY),'-o')

plot(x,-V(:,TRY),'-o')
 hold on 


 xc=linspace(xint,xend,20*N);
SQK=((TRY+0.5).*pi)./2;
DBK=SQK.^2;
 
 COMP=1*sin(xc.*sqrt(Omegalist(TRY,1)));
 COMP2=1*cos(xc.*sqrt(Omegalist(TRY,1)));
%  COMP=1*sin(xc.*SQK);
%  COMP2=1*cos(xc.*SQK);
 plot(xc,COMP,'LineWidth', 1)
 hold on
 title('TRY=',TRY)
 legend
%  plot(x,imag(V(:,TRY)))
%  hold on
%  plot(real(Omegalist),imag(Omegalist),'o')
%  hold on

%  for i = 1:length(Omegalist)
%     text(realOme(i), imagOme(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
% end
% 
% hold off;           
            

%%
%det(A-lambda*B)=0 verification
clc
for i=1:N
% TRYY=499;
TRYY=i;
K_single=Omegalist(TRYY,1);
Res=LHS*V(:,TRYY)-K_single.*RHS*V(:,TRYY);
% Bstar=(LHS-K_single.*RHS)*V(:,TRYY);
% Astar=LHS-K_single.*RHS;
% AS=rref(Astar);
RES(i,1)=max(Res);
% det(Astar)
% det(AS)
% max(Res)
end
plot(LST,RES,'o')
hold on
%%





            

        