#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from scipy.special import jv as bessel
from scipy.integrate import quad,quadrature,fixed_quad
from scipy.interpolate import splrep,splev
from tools import save,load,checkdir
from StrongCoupling import StrongCoupling
from CollinearPDF import Cteq6PDF,ToyPDF
from hankel.hankel import HankelTransform

class TMDPDF(object):

  def __init__(self,CPDF):

    self.D={}
    self.CPDF=CPDF 

    #self.SC=StrongCoupling('exact')
    self.SC=StrongCoupling('one loop')
    self.setup()
    self.set_NP_params()

  def setup(self):
    D=self.D
    D['CF']=4.0/3.0
    D['TF']=1.0/2.0
    D['gamma']=0.577215664901532

    D['AA']=9.0/4/np.pi
    D['LamQCD2']=0.2123**2
    D['A']=4*np.pi/9.0


  def set_NP_params(self,Set=1):
    D=self.D

    Q0=1.3
    b0=0.86*Q0

    if Set==1:
      D['bT2max']=0.5**2
      D['Q02']=Q0**2
      D['C3']=b0
      D['zetaF0']=4*D['Q02']
      D['a1']= 0.21
      D['a2']= 0.68
      D['a3']=-0.6*0.21

    elif Set==2:
      D['bT2max']=1.5**2
      D['Q02']=Q0**2
      D['C3']=b0
      D['zetaF0']=4*D['Q02']
      D['a1']= 0.201
      D['a2']= 0.184
      D['a3']=-0.026

    D['NP model']=lambda x,mu2,bT2: np.exp(-0.5*(D['a1']+D['a2']*\
                  np.log(mu2**0.5/3.2) + D['a3']*np.log(100*x))*bT2)

  def get_bstar2(self,bT2):
    return bT2/(1+bT2/self.D['bT2max'])

  def get_mub2(self,bT2):
    C12=4*np.exp(-2*self.D['gamma'])
    bstar2=self.get_bstar2(bT2)
    return C12/bstar2

  # part A

  def integrand4A(self,r):
    D=self.D
    jz=1-D['x']
    z=D['x'] + r * jz
    return D['DI'](z)*D['pdf'](D['x']) + D['DII'](z)*D['pdf'](D['x']/z)

  def get_A(self,x,mub2,zetaF,bT2,flav):

    # Eqs. A10 & A11 of PRD83 114042

    D=self.D

    factor1 = 0.5*np.log(4/(mub2*bT2))-D['gamma']
    factor2 =-0.5*(np.log(bT2*mub2)-2*(np.log(2)-D['gamma']))**2\
             -(np.log(bT2*mub2) - 2*(np.log(2)-D['gamma'])) * np.log(zetaF/mub2)

    if flav.startswith('g'):

      C0 = 0
      C1 = lambda z: 0
      C2 = lambda z: D['TF']*(2*(1-2*x*(1-x))*factor1 + 2*x*(1-x))
      C3 = factor2

    else: 

      C0 = 1
      C1 = lambda z: D['CF']*2*factor1*2
      C2 = lambda z: D['CF']*(2*factor1*(-1-z) + 1-z)
      C3 = factor2

    # get alpha strong
    alphaS=self.SC.get_alphaS(mub2)

    D['DI']  = lambda z: C0 +alphaS/(2*np.pi)*(-(1-x)*C1(z)/z/(1-z) + C1(1)*np.log(1-x) + C3)
    D['DII'] = lambda z: alphaS/(2*np.pi)*(1-x)*(C1(z)/z/(1-z) + C2(z)/z)
    D['pdf'] = lambda x: self.CPDF.get_CPDF(flav,x,D['muF2'])
    D['x']   = x

    func=np.vectorize(self.integrand4A)
    return fixed_quad(func,0,1,n=40)[0]
    #return quadrature(self.integrand4A,0,1)[0]
    #return quad(self.integrand4A,0,1)[0]

  # part B

  def get_B1(self,mub2,bT2,zetaF):
    D=self.D
    alphaS=self.SC.get_alphaS(mub2)
    Kt=-alphaS*D['CF']/np.pi*(np.log(mub2*bT2)-np.log(4)+2*D['gamma'])
    return 0.5*np.log(zetaF/mub2)*Kt

  def integrand4B2(self,mup2):
    D=self.D
    alphaS=self.SC.get_alphaS(mup2)
    gammaF = alphaS*D['CF']/np.pi*(1.5)
    gammaK = 2*alphaS*D['CF']/np.pi
    return 0.5*(gammaF-0.5*np.log(D['zetaF']/mup2)*gammaK)/mup2

  def get_B2(self,mub2,mu2,zetaF):
    #func=np.vectorize(self.integrand4B2)
    #return fixed_quad(func,mub2,mu2,  n=40)[0]
    #return quadrature(self.integrand4B,mub2,mu2)[0]
    #return quad(self.integrand4B2,mub2,mu2)[0]

    D=self.D
    # NEW 130415
    return 2*D['A']/np.pi*(\
      np.log(np.log(mu2/D['LamQCD2'])/np.log(mub2/D['LamQCD2']))\
     -2.0/3.0*np.log(zetaF/D['LamQCD2'])*\
      np.log(np.log(mu2/D['LamQCD2'])/np.log(mub2/D['LamQCD2']))\
      +2.0/3.0*np.log(mu2/mub2)) 

    #return (1/D['AA']/np.pi)*(\
    #         np.log(np.log(mu2/D['LamQCD2'])/np.log(mub2/D['LamQCD2']))\
    #      - (4.0/3.0)*0.5*np.log(mu2/D['LamQCD2'])\
    #      *(np.log(np.log(mu2/D['LamQCD2'])/np.log(mub2/D['LamQCD2'])))\
    #      + (4.0/3.0)*0.5*np.log(mu2/mub2))  

  def get_B(self,mub2,mu2,zetaF,bT2):
    B1=self.get_B1(mub2,bT2,zetaF)
    B2=self.get_B2(mub2,mu2,zetaF)
    return np.exp(B1+B2) 

  # part C

  def get_C(self,x,zetaF,bT2,flav):
    D=self.D
    return D['NP model'](x,zetaF,bT2)

  # combined 

  def get_PDF_bT_space(self,x,bT2,mu2,zetaF,flav):
    D=self.D
    mub2=self.get_mub2(bT2)
    bT2star=self.get_bstar2(bT2)
    D['muF2']=D['C3']**2/bT2star
    A=self.get_A(x,mub2,zetaF,bT2star,flav)
    B=self.get_B(mub2,mu2,zetaF,bT2star)
    C=self.get_C(x,zetaF,bT2,flav)
    return A*B*C

  def get_PDF(self,x,mu2,zetaF,qT,flav):
    integrand=lambda bT: bT*bessel(0,bT*qT)*\
                        self.get_PDF_bT_space(x,bT**2,mu2,zetaF,flav)
    tgral=quad(integrand,1e-4,20)[0]
    return tgral

  def get_PDF_FFT(self,x,mu2,zetaF,qT,flav):
    pdf=lambda bT: self.get_PDF_bT_space(x,bT**2,mu2,zetaF,flav)
    f=lambda x: x*pdf(x/qT)
    f=np.vectorize(f)
    h = HankelTransform(nu=0,N=120,h=0.003)  
    return h.transform(f,ret_err=False)[0]/qT**2/(2*np.pi) 

if __name__== "__main__":

  CPDF=Cteq6PDF('./')

  tmd=TMDPDF(CPDF)

  x=0.5
  mu2=91.0**2
  zetaF=mu2
  bT2=2
  qT=5.0
  flav='u'
  tmd.D['zetaF']=mu2

  #print tmd.get_PDF_bT_space(x,1.0,mu2,zetaF,'u') 
  for i in range(100):
    #print tmd.get_PDF(x,mu2,zetaF,qT,flav)
    print tmd.get_PDF_FFT(x,mu2,zetaF,qT,flav)

  #ax=py.subplot(111)
  #BT=np.linspace(1e-3,1.3,100)

  #tmd.set_NP_params(Set=1)
  #BTPDF=[bT*tmd.get_PDF_bT_space(x,bT**2,mu2,zetaF,'u') for bT in BT]
  #ax.plot(BT,BTPDF,label='bTmax=1.3')

  #tmd.set_NP_params(Set=2)
  #BTPDF=[bT*tmd.get_PDF_bT_space(x,bT**2,mu2,zetaF,'u') for bT in BT]
  #ax.plot(BT,BTPDF,label='bTmax=0.5')


  #ax.set_ylim(-0.8,0.1)
  #py.savefig('plot.pdf')

  #tmd.profiler.print_stats()










