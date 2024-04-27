clc
clear

N=40;
x=linspace(1,10,N);
h=(10-1)/(N-1);

K1=1/h^2;
K2=1-2/h^2;
K3=1/h^2;

R1=1/h;
R2=-1/h;
R3=0;

LHS=zeros(N);
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
[V,D]=eigs(LHS,N,'smallestabs');      
Omegalist=diag(D);      

TRY=3;
 plot(x,V(:,TRY),'-o')
 hold on
            
            
            
            
            
            

        