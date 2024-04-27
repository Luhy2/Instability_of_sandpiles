clc
clear

Am1 = 200;
Am2 = 400;
Nspace2 = 151;

A = linspace(Am1,Am2,Nspace2);
h = (Am2-Am1)/(Nspace2-1);

lam1 = -1/(2*h);
lam2 = 1;
lam3 = 1/(h^2);
lam4 = -2/(h^2);

MATSIZE = 2*Nspace2;
LLHS = zeros(MATSIZE);
RRHS = eye(MATSIZE);
lkg=1;
for i=1:MATSIZE
    for j=1:MATSIZE
        if lkg==3
            lkg=1;
        end
        if lkg==1
            if i==1
                LLHS(i,i)=1;
                RRHS(i,i)=0;
            elseif i==MATSIZE-1
                %LLHS(i,i-2)=-1/h;
                %LLHS(i,i)=1/h;
                %LLHS(i,i+1)=1;
                LLHS(i,i)=3/(2*h);
                LLHS(i,i+1)=1;
                LLHS(i,i-2)=-2/h;
                LLHS(i,i-4)=1/(2*h);
            else
                LLHS(i,i-2)=-1/(2*h);
                LLHS(i,i+1)=1;
                LLHS(i,i+2)=1/(2*h);
            end
        elseif lkg==2
            if i==2
                LLHS(i,i)=1;
                RRHS(i,i)=0;
            elseif i==MATSIZE
                LLHS(i,i)=1;
                RRHS(i,i)=0;
            else
                LLHS(i,i-2)=1/h^2;
                LLHS(i,i-1)=1;
                LLHS(i,i)=-2/(h^2);
                LLHS(i,i+2)=1/(h^2);
                
            end
        end
        
    end
    lkg=lkg+1;
end



[V,D]=eig(LLHS,RRHS);
D=D./1i;

Omegalist = diag(D);

RL = real(Omegalist);
IM = imag(Omegalist);
%%
plot(RL,IM,'o')
hold on
for i = 1:1:length(Omegalist)
    text(RL(i), IM(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end
%%
figure()
TRY=53;
for i=1:MATSIZE/2
    AlphaH(i,1)=V((i-1)*2+1,TRY);
    AlphaU(i,1)=V((i-1)*2+2,TRY);
end

subplot(2,1,1)
plot(A(1:end),real(AlphaH),'-o')
hold on
plot(A(1:end),imag(AlphaH))
hold on

subplot(2,1,2)
plot(A(1:end),real(AlphaU),'-o')
hold on
plot(A(1:end),imag(AlphaU))
hold on








