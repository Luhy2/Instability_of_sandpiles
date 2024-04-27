clc
clear

F=1.02;
R=8.45;
Gam=0.05877;

Nspace=600;
Aint=100;
Aend=1000;

A=linspace(Aint,Aend,Nspace);
h=(Aend-Aint)/(Nspace-1);

L1=-1/(2*h);
L2=-1/(2*h);
L3=1/(2*h);
L4=1/(2*h);
L5=-1/(2*(F^2)*h);
L6=-(3/2)*(Gam/(F^2));
L7=2/(R*(h^2))+Gam/(F^2);
L8=1/(2*(F^2)*h);
L9=-1/(2*h)-1/(R*(h^2));
L10=1/(2*h)-1/(R*(h^2));

M1=(3/(2*h))-1i*0.505;
M2=-2/h;
M3=1/(2*h);

LeftSub=[1,0;0,1];
Submat=[L1,L2,0,0,L3,L4;L5,L9,L6,L7,L8,L10];
RightSub=[1,0;0,1];

MATSIZE=2*Nspace;

LHS=zeros(MATSIZE);
RHS=zeros(MATSIZE);

Switch=1;
for i=1:MATSIZE
    for j=1:MATSIZE
        if Switch==3
            Switch=1;
        end
        if i>=3 && i<=MATSIZE-2
            if i==j
                RHS(i,j)=1i;
                if Switch==1
                    LHS(i:i+1,j-2:j+3)=Submat;
                end
            end
        elseif i>MATSIZE-2
            if i==MATSIZE-1
                LHS(i,MATSIZE-5)=M3;
                LHS(i,MATSIZE-3)=M2;
                LHS(i,MATSIZE-1)=M1;
            end
            if i==MATSIZE
                LHS(i,MATSIZE-4)=M3;
                LHS(i,MATSIZE-2)=M2;
                LHS(i,MATSIZE)=M1;
            end
        else
            if i==j
                LHS(i,j)=1;
            end
        end
    end
    Switch=Switch+1;
end

% [V,D]=eigs(LHS,RHS,10,'sm');
[V,D]=eig(LHS,RHS);

Omegalist=diag(D);
%%
% figure()
realOme=real(Omegalist);
imagOme=imag(Omegalist);
plot(realOme,imagOme,'o')
hold on

for i = 1:1:length(Omegalist)
    text(realOme(i), imagOme(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end

%%
TRY=576;

for i=1:MATSIZE/2
    AlphaH(i,1)=V((i-1)*2+1,TRY);
    AlphaU(i,1)=V((i-1)*2+2,TRY);
end
figure()
subplot(2,1,1)
plot(A,real(AlphaH),'-o')
hold on
plot(A,imag(AlphaH))
hold on

subplot(2,1,2)
plot(A,real(AlphaU),'-o')
hold on
plot(A,imag(AlphaU))
hold on
















