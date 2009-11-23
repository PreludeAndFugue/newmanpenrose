import unittest
from sympy import conjugate, S, re, im
from newmanpenrose.equations import *

x = Symbol('x', real = True)
z = Symbol('z')

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
            for constant in (1, 0, -5, S.Zero, S.One):
                self.assertEqual(commutator(constant), 0)

    def test_conjugates(self):
        """com1 is real and com4 is imaginary."""
        print re(com1(x))
        for scalar in (x, z):
            self.assertEqual((conjugate(com1(scalar))
                              - com1(conjugate(scalar))).expand(), 0)
            self.assertEqual((conjugate(com4(scalar))
                              + com4(conjugate(scalar))).expand(), 0)
            #self.assertEqual((com1(scalar) - re(com1(scalar))).expand(), 0)

if __name__ == '__main__':
    unittest.main()