clc
clear

Aint = 0;
Aend = 2*pi;

N = 251;

A = linspace(Aint,Aend,N);
h = (Aend-Aint)/(N-1);

MATSIZE = 2*N;
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
                LLHS(i,i)=1;
                RRHS(i,i)=0;
            else
                LLHS(i,i-2)=1/(h^2);
                LLHS(i,i)=-2/(h^2);
                LLHS(i,i+2)=1/(h^2);
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
%                 LLHS(i,i-1)=0;
                LLHS(i,i)=-2/(h^2);
                LLHS(i,i+2)=1/(h^2);
                
            end
        end
        
    end
    lkg=lkg+1;
end

% 三对角预处理子
P1 = LLHS-tril(LLHS,-2)-triu(LLHS,2); 
% 下三角预处理子
P2 = tril(LLHS,0);


[V,D]=eig(LLHS,RRHS);
% [V,D]=eig(P2\LLHS,P2\RRHS);

D=-D;

Omegalist = diag(D);

RL = real(Omegalist);
IM = imag(Omegalist);
%%
plot(RL,IM,'o')
hold on
xlabel('\omega_R')
ylabel('\omega_I')
% for i = 1:1:length(Omegalist)
%     text(RL(i), IM(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
% end
%%
clc
TRY=296;

compN = 5000;
compx = linspace(Aint,Aend,compN);
RLtry = RL(TRY,1);
compU = sin(compx.*sqrt(RLtry));

 figure()

for i=1:MATSIZE/2
    AlphaH(i,1)=V((i-1)*2+1,TRY);
    AlphaU(i,1)=V((i-1)*2+2,TRY);
end

subplot(2,1,1)
plot(A(1:end),real(AlphaH),'-o')
hold on
plot(A(1:end),imag(AlphaH))
hold on
xlabel('A')
ylabel('H')

subplot(2,1,2)
plot(A(1:end),real(AlphaU),'-o')
hold on
plot(A(1:end),imag(AlphaU))
hold on
% plot(compx,-compU)
xlabel('A')
ylabel('U')





