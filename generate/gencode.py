#!/usr/bin/env python

from sympy import Symbol, diff, exp, simplify, limit, pi, factorial, ratsimp
from sympy.utilities.codegen import codegen

def recursion(m,n,base):
    '''
        Compute expression for rho_m*rho_n starting from the expression `base'
        valid for rho_0*rho_0.
    '''
    if m==0 and n==0: return base
    if m>0:
        new = base - alpha/(m+2)*diff(base,alpha)
        return recursion(m-1,n,new)
    if n>0:
        new = base - beta/(n+2)*diff(base,beta)
        return recursion(m,n-1,new)

def recursiontaylor(m,n,base):
    '''
        Compute expression for rho_m*rho_n starting from the expression `base'
        valid for rho_0*rho_0.
    '''
    if m==0 and n==0: return base
    if m>0:
        new = base - alpha/(m+2)*(diff(base,alpha)-diff(base,epsilon))
        return recursiontaylor(m-1,n,new)
    if n>0:
        new = base - (alpha+epsilon)/(n+2)*diff(base,epsilon)
        return recursiontaylor(m,n-1,new)

if __name__=='__main__':
    mmax = 3
    taylor_thresh = 1e-1
    # Define symbols
    alpha = Symbol("alpha")
    beta  = Symbol("beta")
    a = Symbol("a")
    b  = Symbol("b")
    epsilon  = Symbol("epsilon")
    R = Symbol("R")
    # Coulomb integral for rho_0-rho_0 TODO derive this using SymPy
    a0 = beta**4*(beta**2-3*alpha**2)/(beta**2-alpha**2)**3
    a1 = beta**4*alpha/(beta**2-alpha**2)**2/2
    b0 = alpha**4*(alpha**2-3*beta**2)/(alpha**2-beta**2)**3
    b1 = alpha**4*beta/(alpha**2-beta**2)**2/2
    coulomb = 1/R*(1-(a0+a1*R)*exp(-alpha*R) - (b0+b1*R)*exp(-beta*R) )
    # Overlap integral for rho_0-rho_0 TODO derive this using SymPy
    a0 = 4*alpha**4*beta**4/(alpha**2-beta**2)**3
    a1 = alpha**3*beta**4/(alpha**2-beta**2)**2
    b0 = 4*beta**4*alpha**4/(beta**2-alpha**2)**3
    b1 = beta**3*alpha**4/(beta**2-alpha**2)**2
    overlap = 1/(8*R*pi)*( (a0+a1*R)*exp(-alpha*R) + (b0+b1*R)*exp(-beta*R) )
    # Collect all routines in a list
    functions = []
    baseints = [('cou',coulomb),('olp',overlap)]
    for (name,base) in baseints:
        for m in xrange(mmax+1):
            for n in xrange(m+1):
                print "Computing %s integral for m=%d and n=%d" % (name, m, n)
                formula = recursion(m,n,base)
                formula = (formula.subs([(alpha,1/a),(beta,1/b)]))
                functions.append(("%s_%d_%d"%(name,m,n), formula))
                # Analytical expressions are singular for alpha=beta.
                # When alpha \approx beta, we therefore compute the expression
                # by writing beta=alpha+epsilon and writing a Taylor expansion
                # in epsilon.
                if m==0 and n==0:
                    # The Taylor expansion in epsilon is computed only for the
                    # base integral...
                    form = base.subs(beta,alpha+epsilon)
                    series_exp = 1
                    for itaylor in xrange(1,4*mmax+2):
                        series_exp += (-epsilon*R)**itaylor/factorial(itaylor)
                    form = form.subs(exp(-R*(alpha+epsilon)),exp(-R*alpha)*series_exp)
                    form = ratsimp(form)
                    taylor_expansion = simplify(form.subs(epsilon,0))
                    for itaylor in xrange(4*mmax+1):
                        form = ratsimp(diff(form,epsilon))
                        taylor_expansion += simplify(form.subs(epsilon, 0)*epsilon**(itaylor+1)/factorial(itaylor+1))
                    #print taylor_expansion
                # ...other Taylor expansions are derived using a modified recursion scheme.
                formula_taylor = simplify(recursiontaylor(m,n,taylor_expansion))
                functions.append(("%s_taylor_%d_%d"%(name,m,n), formula_taylor.subs([(alpha,1/a),(epsilon,b)])))
    # Write routines in Fortran source code
    print "Writing source code..."
    codegen(functions, "F95", "../src/libslater",to_files=True,
         argument_sequence=(a,b,R), project='libslater')
    # Write python wrapper
    with open('../src/libslater.py','w') as f:
        f.write('''#!/usr/bin/env python\n\n''')
        f.write("'''\nCompute integrals (coulomb and overlap) between two Slater densities defined by\n")
        f.write("    rho_1 = 1/(4 pi (m+2)! a**(m+3)) r^m  exp(-r/a)\n")
        f.write("    rho_2 = 1/(4 pi (n+2)! b**(n+3)) r^n  exp(-r/b)\n")
        f.write("'''\n\n")
        f.write('''import numpy as np\n\nfrom ext import *\n\n''')
        f.write('''__all__ = [%s,'mmax']\n\nmmax = %d\ntaylor_thresh=%f\n\n'''
                % (','.join(["'slater_%s'"%name for (name,base) in baseints]), mmax,taylor_thresh))
        for (name,base) in baseints:
            f.write('''def slater_%s(a,b,m,n,R):\n'''%name)
            f.write('''    if n>m: m,n,a,b = n,m,b,a\n    if m>mmax: raise NotImplementedError, "Regenerate code with mmax>=%%d"%m\n''')
            for m in xrange(mmax+1):
                for n in xrange(m+1):
                    f.write('''    if m==%d and n==%d:\n        if np.abs(1/a-1/b)<taylor_thresh: return %s_taylor_%d_%d(a,1/b-1/a,R)\n        else: return %s_%d_%d(a,b,R)\n'''%(m,n,name,m,n,name,m,n))
            f.write('\n\n')
