#!/usr/bin/env python
'''
Tests for Slater coulomb and overlap expressions
'''

import os
import numpy as np
from scipy.misc import factorial

from horton.grid import AtomicGrid, CubicSpline, solve_poisson_becke
from horton.units import angstrom, kjmol


class SlaterMultipole(object):
    """
    Represents a Slater unit multipole density:
        rho(r) = sqrt((2l+1)/4pi) * Y_lm(omega) r**n-1 exp(-r/sigma) / (l+n+1)! sigma ** l+n+2 )
    """
    def __init__(self,n,l,m,sigma, gridspec='exp:1e-5:1e2:100:974',lmax=8,origin=np.zeros((3,))):
        assert l<n
        assert l<=lmax
        self.n = n
        self.l = l
        self.m = m
        self.sigma = sigma
        self.atgrid = AtomicGrid( 1 , 1, np.zeros((3,)), agspec=gridspec, random_rotate=False )
        self.lmax = lmax
        self.origin = origin
        self.spline_index = self.l*self.l + self.m

    def compute_a_coeffs(self):
        """
        The electrostatic potential of this density can be written as
            phi(r) = sqrt(4pi/2l+1) Y_lm(omega) (1/r**l+1 + sum_k=-l-1^k=n-1 a_k r*k exp(-r/sigma))
        This method returns the (n+l+1,) array of a_k coefficients
        """
        a = np.zeros( (self.n+self.l+1,) )
        # Values from boundary conditions
        a[0]  = -1.0
        a[-1] = -(2.0*self.l + 1.0) / factorial(self.l+self.n+1) / np.power( self.sigma, self.l  + self.n )
        # All other values
        ks = np.arange( -self.l, self.l )
        # Discriminate between l=0 and other cases
        if self.l == 0:
            a[1:2*self.l+1] = - (self.n - ks)/ (self.n + 1) / factorial(ks+1) / np.power(self.sigma, ks+1)
        else:
            a[1:2*self.l+1] = - 1.0 / factorial(ks+self.l+1) / np.power(self.sigma, ks+self.l+1)
        if self.l < self.n - 1:
        # Downward recursion
            a[-2] = 2*self.sigma*(self.n)*a[-1]
            for i in xrange(3,self.n - self.l + 1):
                a[-i] = 2*self.sigma*(self.n+2-i)*a[-i+1] - self.sigma**2*(self.n+2-i-self.l)*(self.n+self.l-i+3)*a[-i+2]
        return a

    def get_rho_splines(self):
        """
        Get the density as a list of splines
        """
        # Turn the radial part of the density into a cubic spline that belongs to spherical harmonic Y_lm
        radii = self.atgrid.rgrid.radii
        rho_radial = np.sqrt( (2.0*self.l + 1.0 ) / 4.0 / np.pi) * np.power(radii, self.n-1) * np.exp(-radii/self.sigma ) / factorial(self.l + self.n + 1) / np.power(self.sigma,1*self.l+self.n+2)
        rho_spline = CubicSpline( rho_radial, rtransform = self.atgrid.rgrid.rtransform )
        # Make a list of CubicSplines
        rho_splines = []
        for i in xrange( (self.lmax+1)**2):
            if i==self.spline_index:rho_splines.append(rho_spline)
            else: rho_splines.append( CubicSpline( np.zeros( radii.shape), rtransform = self.atgrid.rgrid.rtransform))
        return rho_splines

    def get_rho_grid(self,rho_splines=None,center=None):
        """
        Get the density on the atomic grid
        """
        if rho_splines is None: rho_splines = self.get_rho_splines()
        if center is None: center = self.origin
        # Prepare grid points that will contain rho
        rho_grid = np.zeros( (self.atgrid.size, ))
        self.atgrid.eval_decomposition( rho_splines, center, rho_grid)
        return rho_grid

    def get_phi_grid(self,phi_splines=None,center=None):
        """
        Get the density on the atomic grid
        """
        if phi_splines is None: phi_splines = self.get_phi_splines_numerical()
        if center is None: center = self.origin
        # Prepare grid points that will contain rho
        phi_grid = np.zeros( (self.atgrid.size, ))
        self.atgrid.eval_decomposition( phi_splines, center, phi_grid)
        return phi_grid

    def get_phi_splines_numerical(self):
        """
        Get the potential as a list of splines by numerically solving the Poisson equation
        """
        rho_splines = self.get_rho_splines()
        phi_splines = solve_poisson_becke(rho_splines)
        return phi_splines

    def get_phi_splines_analytical(self):
        """
        Get the potential as a list of splines from the analytical solution
        """
        # Get the a_k coefficients and corresponding powers
        a_ks = self.compute_a_coeffs()
        ks = np.arange( -self.l-1, self.n )
        assert len(a_ks) == len(ks)
        # Turn the radial part of the potential into a cubic spline that belongs to spherical harmonic Y_lm
        radii = self.atgrid.rgrid.radii
        phi_spline = np.zeros( radii.shape )
        phi_spline += np.power( radii, -self.l-1)
        for i in xrange(len(a_ks)):
            phi_spline += a_ks[i]*np.power(radii,ks[i])*np.exp(-radii/self.sigma)
        phi_spline = CubicSpline( phi_spline*np.sqrt(4.0*np.pi/(2.0*self.l+1.0)) , rtransform = self.atgrid.rgrid.rtransform )
        # Make a list of CubicSplines
        phi_splines = []
        for i in xrange( (self.lmax+1)**2):
            if i==self.spline_index:phi_splines.append(phi_spline)
            else: phi_splines.append( CubicSpline( np.zeros( radii.shape), rtransform = self.atgrid.rgrid.rtransform))
        return phi_splines


def coulomb_numerical(slatera,slaterb,results=None):
    '''
    Calculate coulomb interaction between two slaters by numerical integration.
    If the result is already stored in a data file, the integration can be skipped.
    If the result is not present yet, the integration is performed and the result is stored.
    '''
    pars_dtype = [
    ('l1', int), ('m1', int), ('l2', int),
    ('m2', int), ('a', float), ('b', float),
    ('x', float), ('y', float), ('z', float),('E',float)]
    pars_names = ['l1','m1','l2','m2','a','b','x','y','z']
    results_fn = 'slater_numint_results.dat'

    def load_results():
        results = np.zeros( 1, dtype=pars_dtype)
        if os.path.isfile(results_fn):
            with open(results_fn) as f:
                for line in f:
                    words = line.split()
                    result = np.zeros( 1, dtype = pars_dtype)
                    for i in xrange(len(words)):
                        result[pars_dtype[i][0]] = pars_dtype[i][1](words[i])
                    if result['E']==0.0: continue
                    results = np.append(results, result)
        return results

    def write_results(results):
        with open('slater_numint_results.dat','w') as f:
            for result in results:
                if result['E']==0.0: continue
                for i in xrange(4): f.write("%3d"%result[pars_dtype[i][0]])
                for i in xrange(4,10): f.write("%16.12f"%result[pars_dtype[i][0]])
                f.write("\n")

    if results is None:
        # Load results from previous numerical integrations
        results = load_results()
    # Setup parameters for this calculation
    this = np.zeros( 1, dtype = pars_dtype)
    this['l1']=slatera.l
    this['m1']=slatera.m
    this['l2']=slaterb.l
    this['m2']=slaterb.m
    this['a']=slatera.sigma
    this['b']=slaterb.sigma
    this['x']=slaterb.origin[0] - slatera.origin[0]
    this['y']=slaterb.origin[1] - slatera.origin[1]
    this['z']=slaterb.origin[2] - slatera.origin[2]
    # Check if this is done before?
    if np.any( this[pars_names] == results[pars_names]):
        # Yes! Return result from previous calculation
        return results[np.where(this[pars_names] == results[pars_names])[0][0]]['E'], results
    # No! Perform calculation
    this['E'] = slatera.atgrid.integrate( slatera.get_rho_grid(), slaterb.get_phi_grid(slaterb.get_phi_splines_analytical(),center=slaterb.origin))
    # Add this result to previous results and store them all
    results = np.append(results,this)
    write_results(results)
    return this['E'], results


def test_point_multipole():
    '''
    Test slater dipole code when all slater radii are zero by comparing with
    point multipole code.
    '''
    from point_multipoles import point_multipole_site_interaction
    R = np.array([1.0,2.0,3.0])
    # All possible combinations separately, allow for easier debugging...
    for i in xrange(4):
        qa = np.zeros((4,))
        qa[i] = 4.0
        for j in xrange(4):
            qb = np.zeros((4,))
            qb[j] = 5.0
            pm = point_multipole_site_interaction(qa,qb,R)
            slater = sp_slater_coulomb(qa,np.zeros((4,)),qb,np.zeros((4,)),R)
            assert np.abs(pm-slater) < 1e-5, "%d %d %+6.3f %+6.3f"%(i,j,pm,slater)
    # all together now...
    qa = np.random.rand(4)
    qb = np.random.rand(4)
    pm = point_multipole_site_interaction(qa,qb,R)
    slater = sp_slater_coulomb(qa,np.zeros((4,)),qb,np.zeros((4,)),R)
    assert np.abs( pm - slater) < 1e-5


def test_slater_dipole():
    '''
    Test slater coulomb expression up to dipoles by comparing with numerical integration.
    WARNING Accurate numerical integration is slow.
    '''
    # Set-up configuration
    # We need to test different combinations of radii, because the implementation
    # is different when both radii are close
    radii = [(2.0,3.0),(2.0,2.1),(2.0,2.05),(2.0,2.01),(2.0,2.00005),(2.0,2.0)]
    alpha = 2.0
    beta = 3.0
    Na = 4.0
    Nb = 5.0
    R = np.array([1.0,2.0,3.0])
    results = None

    for alpha,beta in radii:
        # Check all combinations of spherical harmonics separately
        for la in xrange(2):
            for ma in xrange(2*la+1):
                # Make the first Slater
                slatera = SlaterMultipole(la+1,la,ma,alpha)
                # Set the parameters in a format suitable for our C++ implementation
                qa = np.zeros((4,))
                qa[la**2+ma] = Na
                alphas = np.zeros((4,))
                alphas[la**2+ma] = alpha
                for lb in xrange(2):
                    for mb in xrange(2*lb+1):
                        # Make the second Slater
                        slaterb = SlaterMultipole(lb+1,lb,mb,beta,origin=R)
                        # Set the parameters in a format suitable for our C++ implementation
                        qb = np.zeros((4,))
                        qb[lb**2+mb] = Nb
                        betas = np.zeros((4,))
                        betas[lb**2+mb] = beta
                        # Perform a numerical integration of density times potential to get the Coulomb interaction energy
                        E_int, results = coulomb_numerical(slatera,slaterb, results=results)
                        E_int *= Na*Nb
                        # Compare with our analytical implementation
                        E_ana = sp_slater_coulomb(qa,alphas,qb,betas,R)
                        debug_output = "%d %d %d %d %8.6f %8.6f %+7.4f %+7.4f %+4.1e %+1.3e"%(la,ma,lb,mb,alpha,beta,E_int,E_ana,((E_ana-E_int)/E_int),(E_int-E_ana)/kjmol)
                        assert np.abs(E_int - E_ana) < 1e-4*kjmol
                        #print debug_output


if __name__=='__main__':
    alpha = 3.0
    beta = 5.0
    R = 4.0
    slatera = SlaterMultipole(2,0,0,alpha,origin=np.array([0.0,0.0,0.0]))
    slaterb = SlaterMultipole(2,0,0,beta,origin=np.array([0.0,0.0,R]))
    E_int, results = coulomb_numerical(slatera,slaterb, results=None)
    print E_int
