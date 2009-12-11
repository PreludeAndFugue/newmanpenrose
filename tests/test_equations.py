import unittest
from sympy import conjugate, S, re, im, I
from newmanpenrose.equations import *

x = Symbol('x', real = True)
y = I*x
z = Symbol('z')

# Test commutator functions
def com1_test(scalar):
    """The Delta-D commutator."""
    return Delta(D(scalar)) - D(Delta(scalar)) - (g + conjugate(g))*D(scalar) \
           - (e + conjugate(e))*Delta(scalar) + (conjugate(t) + p)*delta(scalar) \
           + (t + conjugate(p))*deltab(scalar)

def com2_test(scalar):
    """The delta-D commutator."""
    return delta(D(scalar)) - D(delta(scalar)) - (conjugate(a) + b - conjugate(p))*D(scalar) \
           - k*Delta(scalar) + (conjugate(r) + e - conjugate(e))*delta(scalar) \
           + s*deltab(scalar)

def com3_test(scalar):
    """The delta-Delta commutator."""
    return delta(Delta(scalar)) - Delta(delta(scalar)) + conjugate(n)*D(scalar) \
           - (t - conjugate(a) - b)*Delta(scalar) - (m - g + conjugate(g))*delta(scalar) \
           - conjugate(l)*deltab(scalar)

def com4_test(scalar):
    """The delta-deltab commutator."""
    return deltab(delta(scalar)) - delta(deltab(scalar)) - (conjugate(m) - m)*D(scalar) \
           - (conjugate(r) - r)*Delta(scalar) - (a - conjugate(b))*delta(scalar) \
           + (conjugate(a) - b)*deltab(scalar)

test_comms = (com1_test, com2_test, com3_test, com4_test)

class TestScalars(unittest.TestCase):

    def test_complex_scalars(self):
        for scalar in np_scalars:
            self.assertEqual(scalar.is_real, False)

    def test_L(self):
        self.assertEqual(L, conjugate(L))
        self.assertEqual(L.is_real, True)

    def test_real_Phi(self):
        reals = (p00, p11, p22)
        for p in reals:
            self.assertEqual(p, conjugate(p))
            self.assertEqual(p.is_real, True)

    def test_Phi01(self):
        conjugates = ((p01, p10), (p02, p20), (p12, p21))
        for p, pc in conjugates:
            self.assertEqual(p, conjugate(pc))
            self.assertEqual(pc, conjugate(p))

    def test_Psi(self):
        for psi in weyl_scalars:
            self.assertEqual(psi.is_real, False)


class TestCommutators(unittest.TestCase):

    def test_constants(self):
        """Commutators when applied to a constant should return 0."""
        for commutator in comms:
            for constant in (1, 0, -5, S.Zero, S.One, I):
                self.assertEqual(commutator(constant), 0)

    def test_conjugates(self):
        """com1 is real and com4 is imaginary."""
        for scalar in (x, z):
            self.assertEqual((conjugate(com1(scalar))
                              - com1(conjugate(scalar))).expand(), 0)
            self.assertEqual((conjugate(com4(scalar))
                              + com4(conjugate(scalar))).expand(), 0)
            #self.assertEqual((com1(scalar) - re(com1(scalar))).expand(), 0)

    def test_compare_to_test_coms(self):
        """Compare the class based commutators to the function versions defined
        in this module."""
        for c_class, c_test in zip(comms, test_comms):
            self.assertEqual((c_class(z) - c_test(z)).expand(), 0)

    def test_set_rhs_to_zero(self):
        """Test the subs method of the commutator objects by setting the rhs to
        zero."""
        zero_scalars = {k: 0, t: 0, s: 0, r: 0,
                        e: 0, g: 0, a: 0, b: 0,
                        p: 0, n: 0, m: 0, l: 0}
        for com in comms:
            com.subs(zero_scalars)
            op1, op2 = com.op1, com.op2
            self.assertEqual(com(z) - op1(op2(z)) + op2(op1(z)), 0)


class TestCollection(unittest.TestCase):

    def setUp(self):
        """Create a simple collection."""
        self.eqns = {1: e01, 2: e02}
        self.c1 = Collection({1: e01, 2: e02})

    def test_dict(self):
        self.assertEqual(self.c1.equations, self.eqns)

    def test_get_item(self):
        self.assertEqual(self.c1[1], e01)
        self.assertEqual(self.c1[2], e02)

    def test_iter(self):
        for i, eqn in enumerate(self.c1):
            self.assertEqual(eqn, self.eqns[i + 1])
        self.assertEqual(list(self.c1), self.eqns.values())

    def test_subs_zero(self):
        """When all the scalars are set to zero, the equations in the collection
        should all be zero."""
        zero_scalars = {k: 0, t: 0, s: 0, r: 0, e: 0, g: 0, a: 0, b: 0, p: 0,
                        n: 0, m: 0, l: 0,
                        L: 0,
                        p00: 0, p01: 0, p02: 0, p10: 0, p20: 0, p12: 0, p21: 0,
                        p11: 0, p22: 0,
                        psi0: 0, psi1: 0, psi2: 0, psi3: 0, psi4: 0}
        self.c1.subs(zero_scalars)
        for eqn in self.c1:
            self.assertEqual(eqn, 0)

class TestCollectionInstances(unittest.TestCase):

    pass


if __name__ == '__main__':
    unittest.main()