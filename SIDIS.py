#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from TMDPDF import TMDPDF
from TMDFF import TMDFF
from tools import tex
from CollinearPDF import Cteq6PDF,ToyPDF
from FF.FF import FragFuncs
from scipy.integrate import quad
from scipy.special import jv as bessel
from scipy.integrate import quad,quadrature,fixed_quad
from hankel.hankel import HankelTransform

class SIDIS(object):

  def __init__(self,CPDF,CFF):

    self.D={}
    self.setup_params() 
    self.TPDF=TMDPDF(CPDF)
    self.TFF=TMDFF(CFF)
    self.setup_hankel()

  def setup_params(self):
    D=self.D
    D['alphaEM']=1/137.035999074
    D['eU2']=(2.0/3.0)**2
    D['eD2']=(-1.0/3.0)**2
    D['GeV**-2 -> nb'] = 0.3894e6
    D['sin2(thetaW)'] = 0.2397
    D['cos2(thetaW)'] = 1-0.2397
    D['8 pi**2 alphaEM2']=8*np.pi**2*D['alphaEM']**2

  def setup_hankel(self,N=500,h=1e-5):
    D=self.D
    D['hankel']=HankelTransform(nu=0,N=N,h=h)  

  def set_bT2max(self,bT2max):
    D=self.D
    self.TPDF.D['bT2max']=bT2max
    self.TFF.D['bT2max']=bT2max

  def get_L_bT_space(self,x,y,z,Q2,mu2,zetaF,zetaD,bT2,charge):
    D=self.D

    prefactor=z**2*y*D['8 pi**2 alphaEM2']/Q2/Q2/y*(1-y+y*y/2)

    mu2=Q2
    zetaF=Q2
    zetaD=Q2

    TPDF= lambda x,flav: self.TPDF.get_PDF_bT_space(x,bT2,mu2,zetaF,flav)
    TFF = lambda z,flav: self.TFF.get_FF_bT_space(z,bT2,mu2,zetaD,flav,charge)

    L=0
    L+=D['eU2']*TPDF(x,'u')*TFF(z,'u')
    L+=D['eD2']*TPDF(x,'d')*TFF(z,'d')
    L+=D['eU2']*TPDF(x,'ub')*TFF(z,'ub')
    L+=D['eD2']*TPDF(x,'db')*TFF(z,'db')

    return prefactor*L*D['GeV**-2 -> nb']

  def get_L(self,x,y,z,Q2,mu2,zetaF,zetaD,qT,charge):
    integrand=lambda bT: bT*bessel(0,bT*qT)*self.get_L_bT_space(x,y,z,Q2,mu2,zetaF,zetaD,bT**2,charge)/qT
    #tgral=quad(integrand,1e-4,20)[0]
    tgral=quad(integrand,0,np.inf)[0]
    return tgral*qT/(2*np.pi)

  def get_L_FFT(self,x,y,z,Q2,mu2,zetaF,zetaD,qT,charge):
    D=self.D
    L=lambda bT: self.get_L_bT_space(x,y,z,Q2,mu2,zetaF,zetaD,bT**2,charge)/qT
    f=lambda z: z*L(z/qT)
    f=np.vectorize(f)
    return D['hankel'].transform(f,ret_err=False)[0]/qT/(2*np.pi) 

if __name__== "__main__":

  #sidis=SIDIS()
  CPDF=Cteq6PDF('./')
  CFF=FragFuncs('FF/tables/PILO.TAB')

  Q2=15.0
  z=0.5
  x=0.15
  s=300.0
  y=Q2/x/s
  qT=1.0
  charge=1

  sidis=SIDIS(CPDF,CFF)
  #print sidis.get_L_bT_space(0.5,0.7,0.3,2.0,2.0,+1)


  def PDFNP(x,mu2,bT2):
    arg = -0.25*(0.19*np.log(Q2/(1.0)))*bT2
    arg+= -1.0*((0.013*(1-x)**3)/x)*bT2                                          
    return np.exp(arg)
  
  def FFNP(z,mu2,bT2):
    arg = -0.25*(0.19*np.log(Q2/(1.0)))*bT2
    arg+= -(0.2)*bT2                                           
    return np.exp(arg)
  
  sidis.TPDF.D['NP model'] = PDFNP
  sidis.TFF.D['NP model'] = FFNP


  #print sidis.get_L(x,y,z,Q2,qT,+1)
  mu2=Q2
  zetaF=Q2
  zetaD=Q2
  bT2=2.0
  print sidis.get_L_bT_space(x,y,z,mu2,zetaF,zetaD,bT2,charge)
  #print sidis.get_L_FFT(x,y,z,Q2,qT,+1,N=1e2,h=1e-3)
 










