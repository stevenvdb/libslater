#!/usr/bin/env python

import numpy as np

from libslater import slater_cou, slater_olp, mmax

from numerical import SlaterMultipole, coulomb_numerical

def run_tests(name='coulomb'):
    gridspec='exp:1e-5:1e2:200:974'
    exponents = [(0.25,0.25),(0.25,0.25+1e-4),(0.25,0.25+1e-3),(0.25,0.25+1e-2),(0.25,0.25+1e-1)]
    distances = [1.0,5.0]
    print "Testing %s integrals" % name
    print "%5s %5s | %8s %8s %6s | %8s     %8s     %8s     %8s " % ("m","n","alpha","beta","R","Numeric","Analytic","Error","Rel. Err")
    print "="*100
    for m in xrange(mmax+1):
        for n in xrange(mmax+1):
            for alpha,beta in exponents:
                for R in distances:
                    slatera = SlaterMultipole(m+1,0,0,alpha,origin=np.array([0.0,0.0,0.0]),gridspec=gridspec)
                    slaterb = SlaterMultipole(n+1,0,0,beta,origin=np.array([0.0,0.0,R]),gridspec=gridspec)
                    if name=='coulomb':
                        E_num = slatera.atgrid.integrate( slatera.get_rho_grid(), slaterb.get_phi_grid(slaterb.get_phi_splines_analytical(),center=slaterb.origin))
                        E_ana = slater_cou(alpha,beta,m,n,R)
                    elif name=='overlap':
                        E_num = slatera.atgrid.integrate( slatera.get_rho_grid(), slaterb.get_rho_grid() )
                        E_ana = slater_olp(alpha,beta,m,n,R)
                    else: raise NotImplementedError
                    print "%5d %5d | %8.6f %8.6f %6.2f | %+8.2e    %+8.2e    %+8.2e    %+8.6f%%" % (m,n,alpha,beta,R,E_num,E_ana,E_ana-E_num,(E_ana/E_num-1.0)*100.0)
    print "="*100

if __name__=='__main__':
    run_tests(name='coulomb')
    run_tests(name='overlap')
