clc
clear

xint = 0;
xm1 = 20;
Ntot = 20;

x = linspace(xint,xm1,Ntot);
h1 = (xm1-xint)/(Ntot-1);

N = Ntot;
h = h1;

LHSr=zeros(2*N);
rHSr=zeros(2*N);

M1=3/(2*h);
M2=-1;
M3=4/(2*h);
M4=1/(2*h);


for i=1:2*N
    if i==2
        LHSr(i,1)=1;
    elseif i==1
        LHSr(i,2)=-M1;
        LHSr(i,4)=M3;
        LHSr(i,6)=-M4;
        rHSr(i,i)=-1;
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
            LHSr(i,i-1)=-1/(2*h);
            LHSr(i,i+3)=1/(2*h);
            rHSr(i,i)=-1;
        end
        if mod(i,2)==0
            %even
            LHSr(i,i)=-1/2;
            LHSr(i,i-1)=1/h;
            LHSr(i,i-2)=-1/2;
            LHSr(i,i-3)=-1/h;
        end
    end
end

[Vr,Dr]=eigs(LHSr,rHSr,2*N,'smallestabs');
OmegalistR=diag(Dr);
%%
% staggered mesh method diagrams
clc
LengthR=length(OmegalistR);
LSTr=linspace(1,LengthR,LengthR);

plot(LSTr,OmegalistR(1:end,1),'o')
hold on
%%
%figure()
        
TRY=4;

Vr1 = [];
Vr2 = [];


odd_indices = 1:2:length(Vr(:,TRY));
Vr1 = Vr(odd_indices,TRY);

even_indices = 2:2:length(Vr(:,TRY));
Vr2 = Vr(even_indices,TRY);

plot(x,Vr1,'-o')
 hold on
 
%  plot(x,Vr2,'o')
%  hold on 