import operator
from sympy import Symbol, Function, conjugate, S, Mul, Add

#D = Function('D')
#Delta = Function ('Delta')
#delta = Function('delta')
#deltab = Function('deltab')

class DerivativeOperator(Function):

    nargs = 1

    @classmethod
    def eval(cls, arg):
        # TODO:  find out what the difference is between is_Number and is_number
        if arg.is_Number or arg.is_number:
            return S.Zero
        if arg.is_negative:
            return -cls(-arg)
        if arg.could_extract_minus_sign():
            return -cls(-arg)
        if arg.is_Add:
            parts = []
            for val in arg.args:
                parts.append(cls(val))
            return Add(*parts)
        if arg.is_Pow:
            # assume exponent is a constant
            x, n = arg.as_base_exp()
            return n*x**(n - 1)*cls(x)
        if arg.is_Mul:
            factors = []
            for i, val in enumerate(arg.args):
                parts = arg.args[:i] + arg.args[i + 1:] + (cls(val),)
                factors.append(Mul(*parts))
            return Add(*factors)


class D(DerivativeOperator):

    def _eval_conjugate(self):
        # for real derivative operator
        return self.func(conjugate(self.args[0]))


class Delta(DerivativeOperator):

    def _eval_conjugate(self):
        # for real derivative operator
        return self.func(conjugate(self.args[0]))


class delta(DerivativeOperator):

    def _eval_conjugate(self):
        return deltab(conjugate(self.args[0]))


class deltab(DerivativeOperator):

    def _eval_conjugate(self):
        return delta(conjugate(self.args[0]))


# Create the commutators from a class
class Commutator(object):
    """The base commutator object.

    op1 and op2 are the derivative operators on the left-hand-side of the
    commutator.

    rhs is a list of the four coefficients of the derivative operators on the
    right-hand-side."""

    operators = (D, Delta, delta, deltab)

    def __init__(self, op1, op2, rhs):
        self.op1 = op1
        self.op2 = op2
        self.rhs = rhs

    def __call__(self, scalar):
        """Return the result of applying the commutator to a scalar."""
        rhs_parts = [coeff*op(scalar) for coeff, op in zip(self.rhs, Commutator.operators)]
        rhs = map(operator.add, rhs_parts)
        return op1(op2(scalar)) - op2(op1(scalar)) - rhs
        
    def subs(self, substitutions):
        """Simplify the commutator by providing relationships for replacing
        sympy Symbols in the coefficients on the right-hand-side."""
        self.rhs[:] = [coeff.subs(substitution) for coeff in self.rhs]


def main():
    x = Symbol('x', real = True)
    y = Symbol('y')
    z = Symbol('z')

    print conjugate(D(x**2))
    print 'D(5) =', D(5)
    print 'D(-x) = ', D(-x)
    print 'D(x + y) = ', D(x + y)
    print D(x - y)
    print D(x**2)
    print D(x*y)
    #print D(x*y*z)
    print
    print conjugate(delta(x))
    print D(Delta(x**2))
    print
    print
    print Delta.nargs

if __name__ == '__main__':
    main()