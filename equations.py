from sympy import Symbol, conjugate
from operators import D, Delta, delta, deltab

# The NP scalars
k = Symbol('kappa', real = False)
t = Symbol('tau', real = False)
s = Symbol('sigma', real = False)
r = Symbol('rho', real = False)

e = Symbol('epsilon', real = False)
g = Symbol('gamma', real = False)
b = Symbol('beta', real = False)
a = Symbol('alpha', real = False)

p = Symbol('pi', real = False)
n = Symbol('nu', real = False)
m = Symbol('mu', real = False)
l = Symbol('lambda', real = False)

# perhaps this could be a dictionary
# np_scalars = {'k': k, 't': t, ...}
np_scalars = (k, t, s, r, e, g, b, a, p, n, m, l)

# Ricci scalar
L = Symbol('Lambda', real = True)

# Trace-free Ricci tensor components
p00 = Symbol('Phi00', real = True)
p01 = Symbol('Phi01', real = False)
p02 = Symbol('Phi02', real = False)
# p10 is complex-conjugate of p01
p10 = conjugate(p01)
p11 = Symbol('Phi11', real=True)
p12 = Symbol('Phi12', real = False)
# p20 is complex-conjugate of p02
p20 = conjugate(p02)
# p21 is complex-conjugate of p12
p21 = conjugate(p12)
p22 = Symbol('Phi22', real=True)

# perhaps a dictionary...?
trace_free_ricci = (p00, p01, p02, p10, p11, p12, p20, p21, p22)

# Weyl tensor components
psi0 = Symbol('Psi0', real = False)
psi1 = Symbol('Psi1', real = False)
psi2 = Symbol('Psi2', real = False)
psi3 = Symbol('Psi3', real = False)
psi4 = Symbol('Psi4', real = False)

# perhaps a dictionary...?
weyl_scalars = (psi0, psi1, psi2, psi3, psi4)

# spin coefficient equations
e01 = D(r) - deltab(k) - (r*r + s*conjugate(s)) - (e + conjugate(e))*r \
      + conjugate(k)*t + k*(3*a + conjugate(b) - p) - p00
      
e02 = D(s) - delta(k) - (r + conjugate(r) + 3*e - conjugate(e))*s \
      + (t - conjugate(p) + conjugate(a) + 3*b)*k - psi0
      
e03 = D(t) - Delta(k) - (t + conjugate(p))*r - (conjugate(t) + p)*s \
      - (e - conjugate(e))*t + (3*g + conjugate(g))*k - psi1 - p01
      
e04 = D(a) - deltab(e) - (r + conjugate(e) - 2*e)*a - b*conjugate(s) \
      + conjugate(b)*e + k*l + conjugate(k)*g - (e + r)*p - p10
      
e05 = D(b) - delta(e) - (a + p)*s - (conjugate(r) - conjugate(e))*b + (m + g)*k \
      + (conjugate(a) - conjugate(p))*e - psi1

e06 = D(g) - Delta(e) - (t + conjugate(p))*a - (conjugate(t) + p)*b \
      + (e + conjugate(e))*g + (g + conjugate(g))*e - t*p + n*k - psi2 + L - p11

e07 = D(l) - deltab(p) - (r - 3*e + conjugate(e))*l - conjugate(s)*m \
      - (p + a - conjugate(b))*p + n*conjugate(k) - p20

e08 = D(m) - delta(p) - (conjugate(r) - e - conjugate(e))*m - s*l \
      - (conjugate(p) - conjugate(a) + b)*p + n*k - psi2 - 2*L

e09 = D(n) - Delta(p) - (p + conjugate(t))*m - (conjugate(p) + t)*l \
      - (g - conjugate(g))*p + (3*e + conjugate(e))*n - psi3 - p21

e10 = Delta(l) - deltab(n) + (m + conjugate(m) + 3*g - conjugate(g))*l \
      - (3*a + conjugate(b) + p - conjugate(t))*n + psi4

e11 = delta(r) - deltab(s) - (conjugate(a) + b)*r + (3*a - conjugate(b))*s \
      - (r - conjugate(r))*t - (m - conjugate(m))*k + psi1 - p01

e12 = delta(a) - deltab(b) - m*r + l*s - a*conjugate(a) - b*conjugate(b) + 2*a*b \
      - (r - conjugate(r))*g - (m - conjugate(m))*e + psi2 - L - p11

e13 = delta(l) - deltab(m) - (r - conjugate(r))*n - (m - conjugate(m))*p \
      - (a + conjugate(b))*m - (conjugate(a) - 3*b)*l + psi3 - p21

e14 = Delta(m) - delta(n) + (m + g + conjugate(g))*m + l*conjugate(l) \
      - conjugate(n)*p - (conjugate(a) + 3*b - t)*n + p22

e15 = Delta(b) - delta(g) - (conjugate(a) + b - t)*g + m*t - s*n + e*conjugate(n) \
      - (g - conjugate(g) - m)*b + a*conjugate(l) + p12

e16 = Delta(s) - delta(t) + (m - 3*g + conjugate(g))*s + conjugate(l)*r \
      + (t + b - conjugate(a))*t - k*conjugate(n) + p02

e17 = Delta(r) - deltab(t) - (g + conjugate(g) - conjugate(m))*r + s*l \
      - (conjugate(b) - a -conjugate(t))*t - n*k + psi2 + 2*L

e18 = Delta(a) - deltab(g) - (r + e)*n + (t + b)*l - (conjugate(g) - conjugate(m))*a \
      - (conjugate(b) - conjugate(t))*g + psi3

# perhaps a dictionary...?
field_equations = (e01, e02, e03, e04, e05, e06, e07, e08, e09, e10, e11, e12,
                   e13, e14, e15, e16, e17, e18)

# Bianchi identities
b01 = D(psi1) - deltab(psi0) - D(p01) + delta(p00) - (p - 4*a)*psi0 \
      - 2*(2*r + e)*psi1 + 3*k*psi2 + (conjugate(k) - 2*conjugate(a) - 2*b)*p00 \
      + 2*(conjugate(r) + e)*p01 + 2*s*p10 - 2*k*p11 - conjugate(k)*p02

b02 = Delta(psi0) - delta(psi1) + D(p02) - delta(p01) - (4*g - m)*psi0 \
      + 2*(2*t + b)*psi1 - 3*s*psi2 + conjugate(l)*p00 - 2*(conjugate(p) - b)*p01 \
      - 2*s*p11 - (conjugate(r) + 2*e - 2*conjugate(e))*p02 + 2*k*p12

b03 = D(psi2) - deltab(psi1) + Delta(p00) - deltab(p01) + 2*D(L) + l*psi0 \
      - 2*(p - a)*psi1 - 3*r*psi2 - (2*g + 2*conjugate(g) - conjugate(m))*p00 \
      + 2*(a + conjugate(t))*p01 + 2*t*p10 - 2*r*p11 - conjugate(s)*p02

b04 = Delta(psi1) - delta(psi2) - Delta(p01) + deltab(p02) - 2*delta(L) - n*psi0 \
      - 2*(g - m)*psi1 + 3*t*psi2 - 2*s*psi3 + conjugate(n)*p22 \
      - 2*(conjugate(m) - g)*p01 - (2*a + conjugate(t) - 2*conjugate(b))*p02 \
      - 2*t*p11 + 2*r*p12

b05 = D(psi3) - deltab(psi2) - D(p21) + delta(p20) - 2*deltab(L) + 2*l*psi1 \
      - 3*p*psi2 - 2*(r - e)*psi3 + k*psi4 - 2*m*p10 + 2*p*p11 \
      + (2*b + conjugate(p) - 2*conjugate(a))*p20 + 2*(conjugate(r) - e)*p21 \
      - conjugate(k)*p22

b06 = 1

b07 = 1

b08 = 1

b09 = 1

# contracted Bianchi identities
b10 = 1

b11 = 1

b12 = 1

# perhaps a dictionary...?
bianchi_identities = (b01, b02, b03, b04, b05, b06, b07, b08, b09)

# commutator equations (operators)
def com1(scalar):
    """The Delta-D commutator."""
    return Delta(D(scalar)) - D(Delta(scalar)) - (g + conjugate(g))*D(scalar) \
           - (e + conjugate(e))*Delta(scalar) + (conjugate(t) + p)*delta(scalar) \
           + (t + conjugate(p))*deltab(scalar)
           
def com2(scalar):
    """The delta-D commutator."""
    return delta(D(scalar)) - D(delta(scalar)) - (conjugate(a) + b - conjugate(p))*D(scalar) \
           - k*Delta(scalar) + (conjugate(r) + e - conjugate(e))*delta(scalar) \
           + s*deltab(scalar)
           
def com3(scalar):
    """The delta-Delta commutator."""
    pass
    
def com4(scalar):
    """The delta-deltab commutator."""
    pass
    
def commutators(scalar):
    """All the commutators applied to a scalar."""
    return com1(scalar), com2(scalar), com3(scalar), com4(scalar)

def main():
    print e01
    print
    print e02
    print p00.is_real

if __name__ == '__main__':
    main()