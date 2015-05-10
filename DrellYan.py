#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from TMDPDF import TMDPDF
from tools import tex

class DrellYan(object):

  def __init__(self,CPDF):
    self.D={}
    self.setup_params() 
    self.TMD=TMDPDF(CPDF)

  def setup_params(self):
    D=self.D
    D['alphaEM']=1/137.035999074
    D['eU2']=(2.0/3.0)**2
    D['eD2']=(-1.0/3.0)**2
    D['GeV**-2 -> nb'] = 0.3894e6
    D['sin2(thetaW)'] = 0.2397
    D['cos2(thetaW)'] = 1-0.2397

  def get_born(self,Q2,S):
    D=self.D
    return 4*np.pi**3*D['alphaEM']/(3*S)

  def get_Z_boson_charge(self,eq2):
    D=self.D
    return ((1-4*eq2**0.5*D['sin2(thetaW)'])**2+1)\
          /16.0/D['sin2(thetaW)']/D['cos2(thetaW)']

  def get_L_bT_space(self,Q,RS,y,bT):

    D=self.D

    born = self.get_born(Q**2,RS**2)
    x1 = Q/RS*np.exp( y)
    x2 = Q/RS*np.exp(-y)
    mu2=Q*Q
    xiF=Q*Q
    pdf = lambda x,flav: self.TMD.get_PDF_bT_space(x,bT**2,mu2,xiF,flav)

    # see Eq.1.3 of NPB 250 119 (CSS paper)
    QUP=self.get_Z_boson_charge(D['eU2'])
    QDO=self.get_Z_boson_charge(D['eD2'])

    L=0
    L+=QUP*(pdf(x1,'u')*pdf(x2,'ub')+pdf(x1,'ub')*pdf(x2,'u'))  
    L+=QDO*(pdf(x1,'d')*pdf(x2,'db')+pdf(x1,'db')*pdf(x2,'d')) 
    #L+=QDO*(pdf(x1,'s')*pdf(x2,'sb')+pdf(x1,'sb')*pdf(x2,'s')) 
    #L+=QUP*(pdf(x1,'c')*pdf(x2,'cb')+pdf(x1,'cb')*pdf(x2,'c')) 
    #L+=QDO*(pdf(x1,'b')*pdf(x2,'bb')+pdf(x1,'bb')*pdf(x2,'b')) 

    return born*L*D['GeV**-2 -> nb']

if __name__== "__main__":

  DY=DrellYan()
  
  Q=91.1876
  RS=1900.0
  y=0

  BT=np.linspace(1e-3,1.4,100)[:1]

  ax=py.subplot(111)


  DY.TMD.set_NP_params(Set=2)
  BTW=[]
  for bT in BT: BTW.append(bT*DY.get_L_bT_space(Q,RS,y,bT))
  ax.plot(BT,BTW,'b--',lw=2,label=tex('b_{max}=%0.2f')%(DY.TMD.D['bT2max']**0.5))

  #DY.TMD.set_NP_params(Set=1)
  #BTW=[]
  #for bT in BT: BTW.append(bT*DY.get_L_bT_space(Q,RS,y,bT))
  #ax.plot(BT,BTW,'m-.',lw=2,label=tex('b_{max}=%0.2f')%(DY.TMD.D['bT2max']**0.5))

  #ax.set_ylabel(tex('b_T\\tilde{W}\;(nb/GeV)'),size=20)
  #ax.set_xlabel(tex('b_T(GeV^{-1})'),size=20)
  #ax.legend()
  #py.savefig('plot.pdf')







