#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from hankel import HankelTransform
from scipy.special import jv as bessel
from scipy.integrate import quad,quadrature,fixed_quad


BT=np.linspace(0,10,100)
W=lambda bT: np.exp(-0.5*bT**2)

qT=2.0

integrand=lambda bT: bT*bessel(0,bT*qT)*W(bT)
tgral=quad(integrand,1e-4,np.inf)[0]
print tgral



f=lambda x: x*W(x/qT)
#h = HankelTransform(nu=0,N=120,h=0.03)  
h = HankelTransform(nu=0,N=120,h=0.003)  
print h.transform(f,ret_err=False)[0]/qT**2

#f = lambda x: 1  #Define the input function f(x)
#h = HankelTransform(nu=0,N=120,h=0.03)  #Create the HankelTransform instance
#print h.transform(f)  






