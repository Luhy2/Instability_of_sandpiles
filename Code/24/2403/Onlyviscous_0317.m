clc
clear


Am1 = 200;
Am2 = 400;
Nspace2 = 201;

A = linspace(Am1,Am2,Nspace2);
h = (Am2-Am1)/(Nspace2-1);

k2 = 0;
sc = deg2rad(29);

Nr = Nspace2-2;
Imat = eye(Nspace2);

lam1 = -1/(2*h)*Imat; %-1/2h
% lam1 = lam1(2:end-1,2:end-1);

lam2 = diag(1./A); %1/A
% lam2 = lam2(2:end-1,2:end-1);

lam3 = diag(1i*k2.*(1./(A.*cos(sc)))); %ik2(1/X2)
% lam3 = lam3(2:end-1,2:end-1);

lam4 = (1/h^2)*Imat;%1/h^2
% lam4 = lam4(2:end-1,2:end-1);

lam5 = (-(2/h^2)*Imat-diag((k2./(A.*cos(sc))).^2)-3*Imat);
% lam5 = lam5(2:end-1,2:end-1);

lam6 = 6*Imat;
% lam6 = lam6(2:end-1,2:end-1);


% init LHS and RHS
MATSIZE=(Nspace2)*3;
LHS=zeros(MATSIZE);
RHS=zeros(MATSIZE);

lg=1;%LHS counting switch times (current location)
lkg=1;%counting every 3 rows

for i=1:MATSIZE
    if lkg==4
        lkg=1;
        lg=lg+1;
    end
    for j=1:MATSIZE
        if lg==1 %left B.C.s
            if lkg==1
                LHS(i,1)=1;
            end
            if lkg==2
                LHS(i,2)=1;
            end
            if lkg==3
                LHS(i,3)=1;
            end
        
        elseif lg==Nspace2 %right B.C.s
            if lkg==1
%                 LHS(i,end-2)=1;
                LHS(i,3*lg-8)=1/(2*h);
                LHS(i,3*lg-7)=1/(2*h);
                LHS(i,3*lg-5)=-4/(2*h);
                LHS(i,3*lg-4)=-4/(2*h);
                LHS(i,3*lg-2)=3/(2*h);
                LHS(i,3*lg-1)=3/(2*h);
                LHS(i,3*lg)=lam3(lg,lg);
                RHS(i,i)=1;
            end
            if lkg==2
                LHS(i,end-1)=1;
            end
            if lkg==3
                LHS(i,end)=1;
            end
        else % interior points
            if lkg==1
                LHS(i,3*lg-5)=lam1(lg,lg);
                LHS(i,3*lg-4)=lam1(lg,lg);
                LHS(i,3*lg-2)=lam2(lg,lg);
                LHS(i,3*lg-1)=lam2(lg,lg);
                LHS(i,3*lg)=lam3(lg,lg);
                LHS(i,3*lg+1)=-lam1(lg,lg);
                LHS(i,3*lg+2)=-lam1(lg,lg);
                RHS(i,i)=1;
            end
            if lkg==2
                LHS(i,3*lg-4)=lam4(lg,lg);
                LHS(i,3*lg-2)=lam6(lg,lg);
                LHS(i,3*lg-1)=lam5(lg,lg);
                LHS(i,3*lg+2)=lam4(lg,lg);
            end
            if lkg==3
                LHS(i,3*lg-3)=lam4(lg,lg);
                LHS(i,3*lg)=lam5(lg,lg);
                LHS(i,3*lg+3)=lam4(lg,lg);
                
            end
        end
    end
    lkg=lkg+1;
end

% RHS(4:end-3,4:end-3) = eye(3*Nr);

[V,D] = eig(LHS,RHS);
D = D./1i;

Omegalist = diag(D);

RL = real(Omegalist);
IM = imag(Omegalist);
%%
plot(RL,IM,'o')
hold on
xlabel('\omega_R')
ylabel('\omega_I')
for i = 1:1:length(Omegalist)
    text(RL(i), IM(i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end


%%
% figure()
TRY=2;
for i=1:MATSIZE/3
    AlphaH(i,1)=V((i-1)*3+1,TRY);
    AlphaU(i,1)=V((i-1)*3+2,TRY);
    AlphaV(i,1)=V((i-1)*3+3,TRY);
end

subplot(3,1,1)
plot(A(1:end),real(AlphaH),'-o')
hold on
plot(A(1:end),imag(AlphaH))
hold on
xlabel('A')
ylabel('Hprime')


subplot(3,1,2)
plot(A(1:end),real(AlphaU),'-o')
hold on
plot(A(1:end),imag(AlphaU))
hold on
xlabel('A')
ylabel('Uprime')

subplot(3,1,3)
plot(A(1:end),real(AlphaV),'-o')
hold on
plot(A(1:end),imag(AlphaV))
hold on
xlabel('A')
ylabel('Vprime')
