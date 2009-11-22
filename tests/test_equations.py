import unittest
from sympy import conjugate
from newmanpenrose.equations import *

class TestScalars(unittest.TestCase):

    def test_complex_scalars(self):
        for scalar in np_scalars:
            self.assertEqual(k.is_real, False)

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

if __name__ == '__main__':
    unittest.main()