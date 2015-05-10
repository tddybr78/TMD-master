#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from scipy.special import jv as bessel
from scipy.integrate import quad,quadrature,fixed_quad
from scipy.interpolate import splrep,splev
from tools import save,load,checkdir

class StrongCoupling(object):

  def __init__(self,order=1,mode='one loop'):
    self.D={}
    self.D['order']=order
    self.D['mode']=mode
    self.setup()
    #self.load_interpolator()

  def setup(self):
    # ATTENTION:  a = alphaS/4pi 
    D=self.D

    if D['mode']=='exact':
      
      # masses and input scale alphaS
      D['aZ'] = 0.118 /(4*np.pi)
      D['mc2']=1.3**2
      D['mb2']=4.5**2
      D['mZ2']=91.187**2
      D['mt2']=172.9**2

      # set beta function
      D['beta']=np.zeros((7,3))
      for Nf in range(3,7): 
        D['beta'][Nf,0]=11.0-2.0/3.0*Nf 
        D['beta'][Nf,1]=102-38.0/3.0*Nf 
        D['beta'][Nf,2]=2857.0/2.0-5033.0/18.0*Nf+325.0/54.0*Nf**2 

      # get matchings 
      D['ab']=self.evolve_a(D['mZ2'],D['aZ'],D['mb2'],5)
      D['at']=self.evolve_a(D['mZ2'],D['aZ'],D['mt2'],5)
      D['ac']=self.evolve_a(D['mb2'],D['ab'],D['mc2'],4)

      # Q2 range for precalc
      D['Q2min']=1.0
      D['Q2max']=1e6
  
    if D['mode']=='one loop':  
      # parameters for approximate alphaS
      D['AA']=9.0/4/np.pi
      D['LamQCD2']=0.2123**2

  def beta_func(self,a,Nf):
    betaf = -self.D['beta'][Nf,0]
    if self.D['order']>=1: betaf+=-a*self.D['beta'][Nf,1]
    if self.D['order']>=2: betaf+=-a*self.D['beta'][Nf,2]
    return betaf*a**2

  def evolve_a(self,Q20,a,Q2,Nf):

    # Runge-Kutta as implemented in pegasus  
    LR = np.log(Q2/Q20)/20
    for k in range(20):
      XK0 = LR * self.beta_func(a,Nf)
      XK1 = LR * self.beta_func(a + 0.5 * XK0,Nf)
      XK2 = LR * self.beta_func(a + 0.5 * XK1,Nf)
      XK3 = LR * self.beta_func(a + XK2,Nf)
      a+= (XK0 + 2.* XK1 + 2.* XK2 + XK3) * 0.166666666666666
    return a

  def get_a(self,Q2):

    D=self.D

    # get boundary conditions
    if D['mt2']<Q2:

      Q20,a0,Nf=D['mt2'],D['at'],6

    elif D['mb2']<Q2 and Q2<D['mt2']: 

      Q20,a0,Nf=D['mb2'],D['ab'],5

    elif D['mc2']<Q2 and Q2<D['mb2']: 

      Q20,a0,Nf=D['mc2'],D['ac'],4

    elif Q2<=D['mc2']:

      Q20,a0,Nf=D['mc2'],D['ac'],3

    a=self.evolve_a(Q20,a0,Q2,Nf)

    return a

  def load_interpolator(self):
    checkdir('data')
    if os.path.isfile('alphaS.interpolator'):
      D['cspline']=load('alphaS.interpolator')
    else:
      self.precalc()

  def precalc(self):
    D=self.D
    Q2=np.linspace(D['Q2min'],D['Q2max'],1000)
    a=[self.get_a(q2) for q2 in Q2]
    D['cspline']=splrep(Q2,a)
    save(D['cspline'],'data/alphaS.interpolator')

  def get_alphaS_one_loop(self,Q2):
    D=self.D
    return 1/(D['AA']*np.log(Q2/D['LamQCD2']))

  def get_alphaS(self,Q2):
    if self.D['mode']=='interp':
      return splev(Q2,self.D['cspline'])*4*np.pi
    elif self.D['mode']=='normal':
      return self.get_a(Q2)*4*np.pi
    elif self.D['mode']=='one loop':
      return self.get_alphaS_one_loop(Q2)
    else:
      print 'ERR in Strong Coupling: wrong mode.'
      print 'Aborting :('


