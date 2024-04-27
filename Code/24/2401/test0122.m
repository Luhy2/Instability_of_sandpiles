clc
clear

Num=40;

x=linspace(0,pi/2,Num);
H=(pi/2-0)/(Num-1);
h=H;
k1=-1/h;
k2=1/h;

MAT=zeros(2*Num);

MAT(1,1)=1;
MAT(2,2)=1;
MAT(2,3)=k1;
gg=1;
for i=3:2*Num-2
    
    if mod (i,2)~=1%odd
        MAT(i,gg)=1/2;
        MAT(i,gg+1)=k1;
        MAT(i,gg+2)=1/2;
        MAT(i,gg+3)=k2;
        kk=1;
    else%even
        MAT(i,gg)=k2;
        MAT(i,gg+1)=1/2;
        MAT(i,gg+2)=k1;
        MAT(i,gg+3)=1/2;
        kk=2;
    end
    if kk==2
        gg=gg+2;
    end
end
MAT(Num*2-1,Num*2-3)=1;
MAT(Num*2-1,Num*2-2)=k1;
MAT(Num*2,Num*2-2)=1/2;


[V,D]=eigs(MAT,10);%return 10 max eig value


for i=1:10
    jj=0;
    for j=1:2*Num

        if mod (j,2)==1
            jj=jj+1;
            VV(jj,i)=V(j,i);
        end
    end
end


for i=1:10
    plot(x,VV(:,i))
    hold on
end




