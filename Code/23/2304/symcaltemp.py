from sympy import *

    #test
    #a=sin(pi)
    #print(a)

    #b=I+1
    #print(b)

## define parameters

# A=u_a^*
# B=\frac{1}{\rho\phi}
# C=\alpha
# D=k_a
# E=k_{\theta}
# F=k_b

# H=\frac{1}{h_2}
# X=\lambda 关键字使用其他字母做替换
A,B,C,D,E,F,G,H = symbols('A,B,C,D,E,F,G,H',real=True)
X = symbols('X')

UP = -(A*D**3*2*I + A*D*F**2*2*I + A*D*E**2*H**2*2*I + 3*A*D**2*H*cos(C) + A*F**2*H*cos(C) + cos(C)*A*H*(-(F**2*1*I + H*sin(C)*F + D**2*1*I + H*cos(C)*D)**2)**(1/2) - A*D*H**2*cos(C)**2*1*I + 2*A*D*F*H*sin(C) - A*F*H**2*cos(C)*sin(C)*1*I)/(2*(F**2 - sin(C)*F*H*1*I + E**2*H**2 - cos(C)*D*H*1*I + D**2))


print(latex(re(UP)))