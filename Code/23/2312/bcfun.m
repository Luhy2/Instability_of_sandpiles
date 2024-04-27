% Boundary conditions

function res =bcfun(ya,yb)
    res = [ya(1);ya(2);ya(3);yb(4);yb(5)];
end
% function res =bcfun(ya,yb)
%     res=[ya(1)];
% end

