from sympy import Function, Basic

D = Function('D')
Delta = Function ('Delta')
delta = Function('delta')
deltab = Function('deltab')

class DerivativeOperator(Function):
    
    nargs = 1
    
    @classmethod
    def eval(cls, arg):
        # if arg is a constant return 0
        # if arg is ...
        pass
        
    def _eval_conjugate(self):
        pass
        
    