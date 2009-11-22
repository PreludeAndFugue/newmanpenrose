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

np_scalars = (k, t, s, r, e, g, b, a, p, n, m, l)

# Ricci scalar
L = Symbol('Lambda', real = True)

# Trace-free Ricci tensor components
p00 = Symbol('Phi00', real = True)
p01 = Symbol('Phi01', real = False)
p02 = Symbol('Phi02', real = False)
# p10 is complex-conjugate of p01
p10 = conjugate(p01, real = False)
p11 = Symbol('Phi11', real=True)
p12 = Symbol('Phi12', real = False)
# p20 is complex-conjugate of p02
p20 = conjugate(p02, real = False)
# p21 is complex-conjugate of p12
p21 = conjugate(p12, real = False)
p22 = Symbol('Phi22', real=True)

trace_free_ricci = (p00, p01, p02, p10, p11, p12, p20, p21, p22)

# Weyl tensor components
psi0 = Symbol('Psi0', real = False)
psi1 = Symbol('Psi1', real = False)
psi2 = Symbol('Psi2', real = False)
psi3 = Symbol('Psi3', real = False)
psi4 = Symbol('Psi4', real = False)

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

e06 = 1

e07 = 1

e08 = 1

e09 = 1

e10 = 1

e11 = 1

e12 = 1

e13 = 1

e14 = 1

e15 = 1

e16 = 1

e17 = 1

e18 = 1
      
field_equations = (e01, e02, e03, e04, e05, e06, e07, e08, e09, e10, e11, e12,
                   e13, e14, e15, e16, e17, e18)

# Bianchi identities
b01 = 1

b02 = 1

b03 = 1

b04 = 1

b05 = 1

b06 = 1

b07 = 1

b08 = 1

b09 = 1

# contracted Bianchi identities
b10 = 1

b11 = 1

b12 = 1

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