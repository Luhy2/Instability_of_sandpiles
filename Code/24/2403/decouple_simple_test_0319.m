clc
clear

w1 = 90;
w2 = 90;

N = 200;
x = linspace(0,20*pi,N);


[D,xx] = cheb0(N-1);

xx = (xx+1)*10*pi;


y2 = x.*cos(sqrt(w2).*x)./(2*sqrt(w2));
y1 = sin(sqrt(w1).*x);

yy1 = sin(sqrt(w1).*xx);
yy2 = xx.*cos(sqrt(w2).*xx)./(2*sqrt(w2));

xxx = [x xx'];
yyy = [y1 yy1'];

[xxxx,idx] = sort(xxx);
yyyy = yyy(idx);

subplot(2,1,1)
plot(x,y1,'-o')
hold on
plot(xx,yy1,'-o')
hold on
subplot(2,1,2)
plot(xxxx,yyyy,'-o')
% plot(x,y2,'-^')
% hold on
% plot(xx,yy2,'-^')
% hold on



function [ D,x ] = cheb0( N )

if N==0, 
    D = 0;
    x=1;
    return,
end
x = -cos(pi * (0:N)/N)';
c = [2; ones(N-1,1);2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X - X';
D = (c*(1./c)')./(dX + (eye(N+1)));          % off-diagonal entries
D = D - diag(sum(D'));                       % diagonal entries
end