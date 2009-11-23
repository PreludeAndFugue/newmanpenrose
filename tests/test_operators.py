import unittest
from sympy import conjugate, S, diff
from newmanpenrose.operators import D, Delta, delta, deltab
from newmanpenrose.equations import *

operators = (D, Delta, delta, deltab)
real_ops = (D, Delta)
complex_ops = (delta, deltab)

n = Symbol('n')
x = Symbol('x', real = True)
z = Symbol('z')

class TestOperators(unittest.TestCase):

    def test_constants(self):
        """The derivative operators applied to constants should result in
        zero."""
        constants = (0, 1, 5, -10, S.Zero)
        for op in operators:
            for const in constants:
                self.assertEqual(op(const), 0)

    def test_conjugates(self):
        for op in real_ops:
            self.assertEqual(conjugate(op(x)), op(x))
            self.assertEqual(conjugate(op(z)), op(conjugate(z)))
        self.assertEqual(conjugate(delta(x)), deltab(x))
        self.assertEqual(conjugate(deltab(x)), delta(x))
        self.assertEqual(conjugate(delta(z)), deltab(conjugate(z)))
        self.assertEqual(conjugate(deltab(z)), delta(conjugate(z)))

    def test_powers(self):
        operands = (x**4, x**n, z**3, z**n, x**(-1), x**(-n), z**(-2), z**(-n))
        for op in operators:
            for operand in operands:
                base, exp = operand.as_base_exp()
                self.assertEqual((op(operand) - diff(operand, base)*op(base)).expand(), 0)
                # same test
                self.assertEqual((op(operand) - exp*base**(exp - 1)*op(base)).expand(), 0)

    def test_addition(self):
        y = x + z
        for op in operators:
            self.assertEqual(op(y), op(x) + op(z))

    def test_mult_constant(self):
        for op in operators:
            for operand in (x, z):
                self.assertEqual(op(4*operand), 4*op(operand))
                self.assertEqual(op(S.One*operand), S.One*op(operand))

    def test_minus(self):
        for op in operators:
            for operand in (x, z):
                self.assertEqual(op(-operand), -op(operand))

if __name__ == '__main__':
    unittest.main()
