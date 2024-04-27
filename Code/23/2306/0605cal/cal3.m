clearvars -except QWEHG
clc
%计算grouped data下的\lambda与k_a，k_{\theta}的关系
%由于无解析表达，代入数值计算

%网格
L = 1;
Num = 20;
A = linspace(-10*L,10*L,Num);
T = linspace(-10*L,10*L,Num);

[Aa,Tt]=meshgrid(A,T);
for k =1:49
    Q=QWEHG(k,1);
    W=QWEHG(k,2);
    E=QWEHG(k,3);
    H=QWEHG(k,4);
    G=QWEHG(k,5);

%
    syms z;
    for i = 1:Num
        for j = 1:Num
            RS1=root(z^2*(2*W + A(i)*Q*3i + 2*E*Q) + z*(- 3*(A(i))^2*Q^2 + E^2*Q^2 + W^2 + A(i)*Q*W*4i + 3*E*Q*W + (A(i))^2*G*H + A(i)*E*Q^2*4i + E^2*G*H*(T(j))^2) - (A(i))^3*Q^3*1i + z^3 - 2*(A(i))^2*Q^2*W + E^2*Q^2*W + A(i)*Q*W^2*1i + E*Q*W^2 + A(i)*E^2*Q^3*1i - 2*(A(i))^2*E*Q^3 + (A(i))^3*G*H*Q*1i + A(i)*E*Q^2*W*3i + (A(i))^2*E*G*H*Q + E^2*G*H*(T(j))^2*W + A(i)*E^2*G*H*Q*(T(j))^2*1i, z, 1);
            RS2=root(z^2*(2*W + A(i)*Q*3i + 2*E*Q) + z*(- 3*(A(i))^2*Q^2 + E^2*Q^2 + W^2 + A(i)*Q*W*4i + 3*E*Q*W + (A(i))^2*G*H + A(i)*E*Q^2*4i + E^2*G*H*(T(j))^2) - (A(i))^3*Q^3*1i + z^3 - 2*(A(i))^2*Q^2*W + E^2*Q^2*W + A(i)*Q*W^2*1i + E*Q*W^2 + A(i)*E^2*Q^3*1i - 2*(A(i))^2*E*Q^3 + (A(i))^3*G*H*Q*1i + A(i)*E*Q^2*W*3i + (A(i))^2*E*G*H*Q + E^2*G*H*(T(j))^2*W + A(i)*E^2*G*H*Q*(T(j))^2*1i, z, 2);
            RS3=root(z^2*(2*W + A(i)*Q*3i + 2*E*Q) + z*(- 3*(A(i))^2*Q^2 + E^2*Q^2 + W^2 + A(i)*Q*W*4i + 3*E*Q*W + (A(i))^2*G*H + A(i)*E*Q^2*4i + E^2*G*H*(T(j))^2) - (A(i))^3*Q^3*1i + z^3 - 2*(A(i))^2*Q^2*W + E^2*Q^2*W + A(i)*Q*W^2*1i + E*Q*W^2 + A(i)*E^2*Q^3*1i - 2*(A(i))^2*E*Q^3 + (A(i))^3*G*H*Q*1i + A(i)*E*Q^2*W*3i + (A(i))^2*E*G*H*Q + E^2*G*H*(T(j))^2*W + A(i)*E^2*G*H*Q*(T(j))^2*1i, z, 3);
            RSS1(i,j)=real(double(RS1));
            RSS2(i,j)=real(double(RS2));
            RSS3(i,j)=real(double(RS3));
        end
    end
    D1=RSS1>0;
    D2=RSS2>0;
    D3=RSS3>0;
    d1=isequal(D1,zeros(Num));
    d2=isequal(D2,zeros(Num));
    d3=isequal(D3,zeros(Num));
    if d1~=1 || d2~=1 || d3~=1
        figure;
        %输出
        for i=1:3
            subplot(2,2,i)
            name=eval(['RSS',num2str(i)]);
            surf(Aa,Tt,name);

            zlim([-inf,inf])
            legend('\lambda at qwq','Location','northwest');

            xlabel('k_a')
            ylabel('k_{\theta}')
            zlabel('\lambda')
        end
        title('at',num2str(k))
    else
        disp('non')
    end
end
