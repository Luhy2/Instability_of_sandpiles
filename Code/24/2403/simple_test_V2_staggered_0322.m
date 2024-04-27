clc
clear
clc
clear


N1=200;
N2=200;%dissipation mesh points -1
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
C3 = 10; % 10

FA=1+C1*(1+tanh(C2*(-1+(2*x/(xint+xm1+C3)))));
% FA=ones(1,Ntot);


% x = linspace(0,861.6,400);
%%
clc
clear
Ntot = 100;
xint = 0;
xm1 = 20;

x = linspace(xint,xm1,Ntot);
h1=(xm1-xint)/(Ntot-1);

hlist = diff(x);


% check mesh and interval list
% plot(x,'o')
% hold on
% plot(hlist,'o')
% figure()
% plot(x,FA,'o')
% hold on
%%

LHSr = zeros(2*Ntot);
rHSr = zeros(2*Ntot);

K1 = -1/h1;
K2 = 1/h1;
K3 = -1;
K4 = 0;

for i = 1:2*Ntot
    if i==1
        LHSr(i,i)=1;
    elseif i==2*Ntot
%         LHSr(i,i-1)=1;

%         LHSr(i,i)=1;

%           LHSr(i,i)=3;
%           LHSr(i,i-2)=-4;
%           LHSr(i,i-4)=1;
        LHSr(i,i) = 1/2;
        LHSr(i,i-2) = 1/2;
    else
        if mod(i,2)==0
            % even
            LHSr(i,i)=K1;
            LHSr(i,i+2)=K2;
            rHSr(i,i+1)=-1;
            
        end
        if mod(i,2)==1
            % odd
            LHSr(i,i)=K2;
            LHSr(i,i+1)=K4;
            LHSr(i,i-1)=K3;
            LHSr(i,i-2)=K1;
        end
    end
end

[Vr,Dr]=eigs(LHSr,rHSr,2*Ntot,'smallestabs');
% [Vr,Dr]=eig(LHSr,rHSr);
OmegalistR=diag(Dr);

%%

% LengthR=length(OmegalistR);
% LSTr=linspace(1,LengthR,LengthR);
realOme=real(OmegalistR);
imagOme=imag(OmegalistR);

plot(realOme,imagOme,'o')
hold on
xlabel('\omega_R')
ylabel('\omega_I')

for i = 1:1:length(OmegalistR)
    text(realOme(i), imagOme(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end
%%
%figure()
clc
MATSIZE=2*Ntot;
TRY=3;
for i=1:(MATSIZE/2)
    phi1(i,1)=Vr((i-1)*2+1,TRY);
    phi2(i,1)=Vr((i-1)*2+2,TRY);
    
end

% figure()
% 
% subplot(2,1,1)
plot(x(1:end),-real(phi1),'-o')
hold on
% legend('grid point','theoretical','reduction','staggered')
legend('grid point','theoretical','staggered')
% plot(x(1:end),imag(phi1))
% hold on



% subplot(2,1,2)
% plot(x(1:end),real(phi2),'-o')
% hold on
% plot(x(1:end),imag(phi2))
% hold on







            

        