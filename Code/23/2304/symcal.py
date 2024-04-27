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
A,B,C,D,E,F,H = symbols('A,B,C,D,E,F,H',real=True)
X = symbols('X')

    #test2
    #c=X+I*D*A+sin(C)*H
    #print(c)

M = Matrix([[X+I*D*A,0,0,B*I*D],[0,X+I*D*A+cos(C)*H*A,0,B*H*I*E],[0,0,X+I*D*A,B*I*F],[I*D+cos(C)*H,I*H*E,I*F+sin(C)*H,0]])

R = M.det()
#print(R)
#print(latex(R))

#此时R需要从字符转回sympy格式
ExprR = sympify(R)
#print(ExprR)


#Re = solve(ExprR,X)
#print(latex(Re))

#R2 = 'A*(-2*I*D**3 - 4*D**2*H*cos(C) - 2*I*D*E**2*H**2 - 2*I*D*F**2 - 2*D*F*H*sin(C) + 2*I*D*H**2*cos(C)**2 - 2*F**2*H*cos(C) + I*F*H**2*sin(2*C))/(2*(D**2 - I*D*H*cos(C) + E**2*H**2 + F**2 - I*F*H*sin(C)))'
#RX2 = sympify(R2)
RX2 = A*(-2*I*D**3 - 4*D**2*H*cos(C) - 2*I*D*E**2*H**2 - 2*I*D*F**2 - 2*D*F*H*sin(C) + 2*I*D*H**2*cos(C)**2 - 2*F**2*H*cos(C) + I*F*H**2*sin(2*C))/(2*(D**2 - I*D*H*cos(C) + E**2*H**2 + F**2 - I*F*H*sin(C)))
#print(RX2)

REX2 = re(RX2)
IMX2 = im(RX2)
#print(latex(REX2))

SrX2 = simplify(REX2)
SiX2 = simplify(IMX2)

#print(latex(SrX2))
#print(latex(SiX2))

print(SiX2)