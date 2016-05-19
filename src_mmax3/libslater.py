#!/usr/bin/env python

'''
Compute integrals (coulomb and overlap) between two Slater densities defined by
    rho_1 = 1/(4 pi (m+2)! a**(m+3)) r^m  exp(-r/a)
    rho_2 = 1/(4 pi (n+2)! b**(n+3)) r^n  exp(-r/b)
'''

import numpy as np

from ext import *

__all__ = ['slater_cou','slater_olp','mmax']

mmax = 3
taylor_thresh=0.100000

def slater_cou(a,b,m,n,R):
    if n>m: m,n,a,b = n,m,b,a
    if m>mmax: raise NotImplementedError, "Regenerate code with mmax>=%%d"%m
    if m==0 and n==0:
        if np.abs(1/a-1/b)<taylor_thresh: return cou_taylor_0_0(a,1/b-1/a,R)
        else: return cou_0_0(a,b,R)
    if m==1 and n==0:
        if np.abs(1/a-1/b)<taylor_thresh: return cou_taylor_1_0(a,1/b-1/a,R)
        else: return cou_1_0(a,b,R)
    if m==1 and n==1:
        if np.abs(1/a-1/b)<taylor_thresh: return cou_taylor_1_1(a,1/b-1/a,R)
        else: return cou_1_1(a,b,R)
    if m==2 and n==0:
        if np.abs(1/a-1/b)<taylor_thresh: return cou_taylor_2_0(a,1/b-1/a,R)
        else: return cou_2_0(a,b,R)
    if m==2 and n==1:
        if np.abs(1/a-1/b)<taylor_thresh: return cou_taylor_2_1(a,1/b-1/a,R)
        else: return cou_2_1(a,b,R)
    if m==2 and n==2:
        if np.abs(1/a-1/b)<taylor_thresh: return cou_taylor_2_2(a,1/b-1/a,R)
        else: return cou_2_2(a,b,R)
    if m==3 and n==0:
        if np.abs(1/a-1/b)<taylor_thresh: return cou_taylor_3_0(a,1/b-1/a,R)
        else: return cou_3_0(a,b,R)
    if m==3 and n==1:
        if np.abs(1/a-1/b)<taylor_thresh: return cou_taylor_3_1(a,1/b-1/a,R)
        else: return cou_3_1(a,b,R)
    if m==3 and n==2:
        if np.abs(1/a-1/b)<taylor_thresh: return cou_taylor_3_2(a,1/b-1/a,R)
        else: return cou_3_2(a,b,R)
    if m==3 and n==3:
        if np.abs(1/a-1/b)<taylor_thresh: return cou_taylor_3_3(a,1/b-1/a,R)
        else: return cou_3_3(a,b,R)


def slater_olp(a,b,m,n,R):
    if n>m: m,n,a,b = n,m,b,a
    if m>mmax: raise NotImplementedError, "Regenerate code with mmax>=%%d"%m
    if m==0 and n==0:
        if np.abs(1/a-1/b)<taylor_thresh: return olp_taylor_0_0(a,1/b-1/a,R)
        else: return olp_0_0(a,b,R)
    if m==1 and n==0:
        if np.abs(1/a-1/b)<taylor_thresh: return olp_taylor_1_0(a,1/b-1/a,R)
        else: return olp_1_0(a,b,R)
    if m==1 and n==1:
        if np.abs(1/a-1/b)<taylor_thresh: return olp_taylor_1_1(a,1/b-1/a,R)
        else: return olp_1_1(a,b,R)
    if m==2 and n==0:
        if np.abs(1/a-1/b)<taylor_thresh: return olp_taylor_2_0(a,1/b-1/a,R)
        else: return olp_2_0(a,b,R)
    if m==2 and n==1:
        if np.abs(1/a-1/b)<taylor_thresh: return olp_taylor_2_1(a,1/b-1/a,R)
        else: return olp_2_1(a,b,R)
    if m==2 and n==2:
        if np.abs(1/a-1/b)<taylor_thresh: return olp_taylor_2_2(a,1/b-1/a,R)
        else: return olp_2_2(a,b,R)
    if m==3 and n==0:
        if np.abs(1/a-1/b)<taylor_thresh: return olp_taylor_3_0(a,1/b-1/a,R)
        else: return olp_3_0(a,b,R)
    if m==3 and n==1:
        if np.abs(1/a-1/b)<taylor_thresh: return olp_taylor_3_1(a,1/b-1/a,R)
        else: return olp_3_1(a,b,R)
    if m==3 and n==2:
        if np.abs(1/a-1/b)<taylor_thresh: return olp_taylor_3_2(a,1/b-1/a,R)
        else: return olp_3_2(a,b,R)
    if m==3 and n==3:
        if np.abs(1/a-1/b)<taylor_thresh: return olp_taylor_3_3(a,1/b-1/a,R)
        else: return olp_3_3(a,b,R)


