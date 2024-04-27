clc
clear

tint=1;
tend=10;
amesh=linspace(tint,tend,100);
solinit=bvpinit(amesh,@guess);

opts = bvpset('RelTol',1e-9,'Stats','on');%指定相对误差容限
sol=bvp5c(@aODE,@bcfun,solinit,opts);

plot(sol.x,sol.y(1,:),'-o')
hold on