clc
clear


N1=200;
N2=200;%dissipation mesh points -1 变疏 90 加密 600
Ntot=N1+N2;

xint=0;
xm1=20;

x1=linspace(xint,xm1,N1);
h1=(xm1-xint)/(N1-1);
xm2=xm1+h1;

alpha=h1/(xm2^2); %increasing rate

x2(1)=xm2;

for i = 2:N2
      x2(i) = x2(i-1)*(1+alpha*x2(i-1));  %变疏
      %x2(i) = x2(i-1)/(1+(alpha)*x2(i-1)); %加密
end

x = [x1 x2]; %变疏
%x = [x1 -x2+2*xm2]; %加密


hlist = diff(x);

C1 = 14; %1 2 14
C2 = 15;
C3 = 400; % 10

% FA=1+C1*(1+tanh(C2*(-1+(2*x/(xint+xm1+C3)))));
FA=ones(1,Ntot);
%%
clc
clear
Ntot = 5;
xint = 0;
xm1 = 20;

x = linspace(xint,xm1,Ntot);
h1=(xm1-xint)/(Ntot-1);
hlist = diff(x);

FA=ones(1,Ntot);
% check mesh and interval list
% figure()
% plot(x,'o')
% hold on
% plot(hlist,'o')
% figure()
% plot(x,FA,'o')
% hold on
%%
L3=zeros(Ntot,1);
L2=zeros(Ntot,1);
L1=zeros(Ntot,1);
for i=2:Ntot-1
    L3(i,1)=2/(hlist(1,i-1)*(hlist(1,i)+hlist(1,i-1))).*FA(1,i);
    L2(i,1)=-2/(hlist(1,i)*hlist(1,i-1)).*FA(1,i);
    L1(i,1)=2/(hlist(1,i)*(hlist(1,i)+hlist(1,i-1))).*FA(1,i);
end

LHS = spdiags([L3 L2 L1],-1:1,Ntot,Ntot);
FLHS = full(LHS);
LLHS = FLHS(2:end-1,2:end-1);

rHS=eye(Ntot);
RHS=-1.*rHS;
RRHS = RHS(2:end-1,2:end-1);

[V,D]=eigs(LLHS,RRHS,Ntot-2,'smallestabs');
%%
Omegalist=diag(D);      
% realOme=real(Omegalist);
% imagOme=imag(Omegalist);
Length=length(Omegalist);
LST=linspace(1,Length,Length);

plot(Omegalist,'o')
hold on
realOme=real(Omegalist);
% imagOme=imag(Omegalist);
% for i = 1:10:length(Omegalist)
%     text(i, realOme(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
% end

%%
clc
% figure()
TRY=1;

% plot(x(1,2:end-1),V(:,TRY),'-o')
plot(x,[0,V(:,TRY)',0],'-o')
 hold on 
% xlim([0,20])
%%
K1=1/h^2;
K2=-2/(h^2);
K3=1/h^2;


R1=3/(2*h);
R2=-2/h;
R3=1/(2*h);

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
%%
% [V,D]=eig(LHS,RHS);
[V,D]=eigs(LHS,RHS,N,'smallestabs');
% [V,D]=eigs(LHS,RHS,N);
Omegalist=diag(D);      
% realOme=real(Omegalist);
% imagOme=imag(Omegalist);
Length=length(Omegalist);
LST=linspace(1,Length,Length);

plot(Omegalist,'o')

%%
% grid point method diagram
clc
%figure()
TRY=11;
%  plot(x,real(V(:,TRY)),'-o')
%  hold on
%  plot(x,-13.23*V(:,TRY),'-o')

plot(x(1,2:end-1),V(:,TRY),'-o')
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



%%
% reduction method
N=Ntot;
h=h1;

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
        LHSr(i,i)=1;
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

[Vr,Dr]=eigs(LHSr,rHSr,2*N,'smallestabs');
OmegalistR=diag(Dr);

%%
% reduction method diagrams
clc
LengthR=length(OmegalistR);
LSTr=linspace(1,LengthR,LengthR);

% plot(LSTr,OmegalistR(1:end,1),'o')
% hold on

%figure()
        
TRY=1;
Vr1 = [];
Vr2 = [];


odd_indices = 1:2:length(Vr(:,TRY));
Vr1 = Vr(odd_indices,TRY);

even_indices = 2:2:length(Vr(:,TRY));
Vr2 = Vr(even_indices,TRY);

plot(x,Vr1,'-o')
 hold on
%  plot(x,VR2,'o')
%  hold on 












            

        