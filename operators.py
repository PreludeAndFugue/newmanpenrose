from sympy import Symbol, Function, conjugate, S

#D = Function('D')
#Delta = Function ('Delta')
#delta = Function('delta')
#deltab = Function('deltab')    

class DerivativeOperator(Function):

    nargs = 1

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            return S.Zero
        if arg.is_negative:
            return -cls(-arg)
        if arg.could_extract_minus_sign():
            return -cls(-arg)
        if arg.is_Add:
            result = 0
            for val in arg.args:
                result += cls(val)
            return result
        if arg.is_Pow:
            x, n = arg.as_base_exp()
            return n*x**(n - 1)*cls(x)
        if arg.is_Mul:
            # TODO: this has to be generalised to arbitrary number of args
            x, y = arg.args
            return x*cls(y) + y*cls(x)


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