from sympy import *

    #test
    #a=sin(pi)
    #print(a)

    #b=I+1
    #print(b)

## define parameters (2)

# A=u_a^*=u_1
# B=\frac{1}{\rho\phi}
# C=\alpha
# D=k_a
# E=k_{\theta}
# F=k_b

# G=u_{\theta}^*=u_2

# H=\frac{1}{h_2}
# X=\lambda 关键字使用其他字母做替换
A,B,C,D,E,F,G,H = symbols('A,B,C,D,E,F,G,H',real=True)
X = symbols('X')

    #test2
    #c=X+I*D*A+sin(C)*H
    #print(c)

M = Matrix([[X+I*D*A+H*I*E*G,-2*H*cos(C)*G,0,B*I*D],[cos(C)*H*G,X+I*D*A+H*I*E*G+cos(C)*H*A,sin(C)*H*G,B*H*I*E],[0,-2*H*cos(C)*G,X+I*D*A+H*I*E*G,B*I*F],[I*D+cos(C)*H,I*H*E,I*F+sin(C)*H,0]])

#print(latex(M))
#print(M)

#R = M.det()

#print(R)
#print(latex(R))

R='-(A*D^3*2*I + (- 8*(F^2 - sin(C)*F*H*1*I + E^2*H^2 - cos(C)*D*H*1*I + D^2)*(2*A^2*D^2*H^2*cos(C)^2 - 2*A^2*D^2*F^2 - 2*E^4*G^2*H^4 - 2*A^2*D^4 + 4*E^2*G^2*H^4*cos(C)^2 + 4*F^2*G^2*H^2*cos(C)^2 + 2*D^2*G^2*H^2*sin(2*C) + 2*E^2*G^2*H^4*sin(2*C) + A^2*D^3*H*cos(C)*4*I - 2*A^2*D^2*E^2*H^2 - 2*D^2*E^2*G^2*H^2 - 2*E^2*F^2*G^2*H^2 - 4*A*D^3*E*G*H + D*E^2*G^2*H^3*cos(C)*4*I - 4*D*F*G^2*H^2*cos(C)^2 + E^2*F*G^2*H^3*cos(C)*4*I + A^2*D*F*H^2*sin(2*C) - 2*D*F*G^2*H^2*sin(2*C) - 4*A*D*E^3*G*H^3 + A^2*D*F^2*H*cos(C)*2*I + A^2*D^2*F*H*sin(C)*2*I + A*D^2*E*G*H^2*cos(C)*8*I + 6*A*D*E*G*H^3*cos(C)^2 + A*E*F^2*G*H^2*cos(C)*2*I + 2*A*D*E*G*H^3*sin(2*C) + A*E*F*G*H^3*sin(2*C) - 4*A*D*E*F^2*G*H + A*D*E*F*G*H^2*cos(C)*4*I + A*D*E*F*G*H^2*sin(C)*2*I) + (A*D^3*4*I + A*D*F^2*4*I + E^3*G*H^3*4*I + A*D*E^2*H^2*4*I + 6*A*D^2*H*cos(C) + 2*A*F^2*H*cos(C) - A*D*H^2*cos(C)^2*2*I - E*G*H^3*cos(C)^2*4*I - A*F*H^2*sin(2*C)*1*I - E*G*H^3*sin(2*C)*2*I + D^2*E*G*H*4*I + E*F^2*G*H*4*I + 4*A*D*F*H*sin(C) + 6*D*E*G*H^2*cos(C) + 4*E*F*G*H^2*cos(C) + 2*E*F*G*H^2*sin(C))^2)^(1/2)/2 + A*D*F^2*2*I + E^3*G*H^3*2*I + A*D*E^2*H^2*2*I + 3*A*D^2*H*cos(C) + A*F^2*H*cos(C) - A*D*H^2*cos(C)^2*1*I - E*G*H^3*cos(C)^2*2*I - (A*F*H^2*sin(2*C)*1*I)/2 - E*G*H^3*sin(2*C)*1*I + D^2*E*G*H*2*I + E*F^2*G*H*2*I + 2*A*D*F*H*sin(C) + 3*D*E*G*H^2*cos(C) + 2*E*F*G*H^2*cos(C) + E*F*G*H^2*sin(C))/(2*(F^2 - sin(C)*F*H*1*I + E^2*H^2 - cos(C)*D*H*1*I + D^2))'
#此时R需要从字符转回sympy格式
ExprR = sympify(R)
print(ExprR)


#Re = solve(ExprR,X)
#print(latex(Re))
#print(Re)

