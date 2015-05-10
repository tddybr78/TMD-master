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
from StrongCoupling import StrongCoupling

class SIDIS_L(object):

  def __init__(self,CPDF,CFF):

    self.D={}
    self.setup_params() 
    self.TPDF=TMDPDF(CPDF)
    self.TFF=TMDFF(CFF)
    self.setup_hankel()
    self.CPDF=CPDF
    self.CFF=CFF

  def setup_params(self):
    D=self.D
    D['alphaEM']=1/137.035999074
    D['alphaEM2']=D['alphaEM']**2
    D['eU2']=(2.0/3.0)**2
    D['eD2']=(-1.0/3.0)**2
    D['GeV**-2 -> nb'] = 0.3894e6
    D['sin2(thetaW)'] = 0.2397
    D['cos2(thetaW)'] = 1-0.2397

  def setup_hankel(self,N=500,h=1e-5):
    D=self.D
    D['hankel']=HankelTransform(nu=0,N=N,h=h)  

  def set_bT2max(self,bT2max):
    D=self.D
    self.TPDF.D['bT2max']=bT2max
    self.TFF.D['bT2max']=bT2max

  def get_FD_bT_space(self,x,y,z,Q2,mu2,zetaF,zetaD,bT2,charge):

    # see Eq.74 of TMD reference

    D=self.D

    TPDF= lambda x,flav: self.TPDF.get_PDF_bT_space(x,bT2,mu2,zetaF,flav)
    TFF = lambda z,flav: self.TFF.get_FF_bT_space(z,bT2,mu2,zetaD,flav,charge)

    FD =D['eU2']*TPDF(x,'u')*TFF(z,'u')
    FD+=D['eD2']*TPDF(x,'d')*TFF(z,'d')
    FD+=D['eU2']*TPDF(x,'ub')*TFF(z,'ub')
    FD+=D['eD2']*TPDF(x,'db')*TFF(z,'db')
    return FD

  def get_FD_qT_space(self,x,y,z,Q2,mu2,zetaF,zetaD,qT,charge,method='FFT'):
    D=self.D

    # hint: 
    #   hankel transform for F*D (see Eq.75 of reference)
    #   FD_qT_space  = \int_0^inf dbT bT J0(qT*bT) F(bT) D(bT)   <- change var: k = qT*bT  
    #                = (1/qT**2) * \int_0^inf dk J0(k) [k F(k/bT) D(k/qT)]

    FD_bT_space=lambda bT: self.get_FD_bT_space(x,y,z,Q2,mu2,zetaF,zetaD,bT**2,charge)
    f=lambda k: k*FD_bT_space(k/qT)

    if method=='FFT':

      f=np.vectorize(f)
      FD_qT_space = D['hankel'].transform(f,ret_err=False)[0] / qT**2 * 2*np.pi

    elif method=='quad':

      integrand=lambda k: bessel(0,k)*f(k)
      FD_qT_space = quad(integrand,1e-4,np.inf)[0] / qT**2 *2*np.pi

    return FD_qT_space

  def get_L(self,x,y,z,Q2,mu2,zetaF,zetaD,qT,charge,ps='dxdQ2dzdqT2',method='FFT'):

    FD_qT_space = self.get_FD_qT_space(x,y,z,Q2,mu2,zetaF,zetaD,qT,charge,method=method)

    D=self.D
    if ps=='dxdydzdqT': 
      prefactor=4*np.pi*z**2*qT*D['alphaEM2']/Q2/y*(1-y+y*y/2) * D['GeV**-2 -> nb']
      return prefactor * FD_qT_space
    elif ps=='dxdQ2dzdqT2': 
      prefactor=4*np.pi*z**2*qT*D['alphaEM2']/Q2/y*(1-y+y*y/2) * D['GeV**-2 -> nb'] * y/Q2/2/qT 
      return prefactor * FD_qT_space
    else: 
      print 'ps not inplemented'
      return None

class SIDIS_FO(object):

  def __init__(self,CPDF,CFF):

    self.D={}
    self.setup_params() 
    self.CPDF=CPDF
    self.CFF=CFF
    self.SC=StrongCoupling('one loop')

  def setup_params(self):
    D=self.D
    D['alphaEM']=1/137.035999074
    D['eU2']=(2.0/3.0)**2
    D['eD2']=(-1.0/3.0)**2
    D['GeV**-2 -> nb'] = 0.3894e6
    D['sin2(thetaW)'] = 0.2397
    D['cos2(thetaW)'] = 1-0.2397
    D['CF']=4.0/3.0

  def get_sum_Xi_A(self,channel,xh,zh,Q2,qT2):
    D=self.D
    if channel=='quark->quark':
      braket= ((Q2/xh/zh)**2+(Q2-qT2)**2)/qT2+6*Q2
      return 2*D['CF']*xh*zh*(braket*D['A1']+4*Q2*D['A2'])
    elif channel=='quark->gluon':
      braket= Q2**2/qT2*(1/(xh*zh)**2-2/(xh*zh)+2) + 2*Q2*(5-1/xh-1/zh)
      return xh*(1-xh)*(braket*D['A1']+8*Q2*D['A2'])
    elif channel=='gluon->quark':
      qT2_=(zh/(1-zh))**2*qT2
      braket= (Q2**2/(xh*(1-zh))**2 + (Q2-qT2_)**2)/qT2_ + 6*Q2
      return 2*D['CF']*xh*(1-zh)*(braket*D['A1']+4*Q2*D['A2'])

  def get_PDF_FF(self,channel,xia,xib,muF2,charge):
    D=self.D
    CPDF=self.CPDF
    CFF=self.CFF
    if channel=='quark->quark':
      out =D['eU2']*CPDF.get_CPDF('u' ,xia,muF2) * CFF.get_FF(xib,muF2,'u' ,charge)
      out+=D['eD2']*CPDF.get_CPDF('d' ,xia,muF2) * CFF.get_FF(xib,muF2,'d' ,charge)
      out+=D['eU2']*CPDF.get_CPDF('ub',xia,muF2) * CFF.get_FF(xib,muF2,'ub',charge)
      out+=D['eD2']*CPDF.get_CPDF('db',xia,muF2) * CFF.get_FF(xib,muF2,'db',charge)
    elif channel=='quark->gluon':
      out =D['eU2']*CPDF.get_CPDF('u' ,xia,muF2) 
      out+=D['eD2']*CPDF.get_CPDF('d' ,xia,muF2) 
      out+=D['eU2']*CPDF.get_CPDF('ub',xia,muF2) 
      out+=D['eD2']*CPDF.get_CPDF('db',xia,muF2) 
      out*=CFF.get_FF(xib,muF2,'g',charge)
    elif channel=='gluon->quark':
      out =D['eU2']*CFF.get_FF(xib,muF2,'u' ,charge) 
      out+=D['eD2']*CFF.get_FF(xib,muF2,'d' ,charge) 
      out+=D['eU2']*CFF.get_FF(xib,muF2,'ub',charge) 
      out+=D['eD2']*CFF.get_FF(xib,muF2,'db',charge) 
      out*=CPDF.get_CPDF('g',xia,muF2)
    return out

  def get_M(self,xia,xib,xh,zh,Q2,muF2,qT2,charge):
    D=self.D
    out=0
    for channel in ['quark->quark','quark->gluon','gluon->quark']:
      out+= self.get_PDF_FF(channel,xia,xib,muF2,charge) * self.get_sum_Xi_A(channel,xh,zh,Q2,qT2)
    out*=xh*zh
    return out

  def get_FO(self,x,y,z,Q2,qT2,muR2,muF2,charge,ps='dxdQ2dzdqT2',method='gauss'):
    D=self.D
    D['A1']=1+(2/y-1)**2
    D['A2']=-2
    w2=qT2/Q2*x*z
    w=w2**0.5
    xia_=lambda xib: w2/(xib-z)+x
    xib_=lambda xia: w2/(xia-x)+z

    integrand_xia=lambda xia: self.get_M(xia,xib_(xia),x/xia,z/xib_(xia),Q2,muF2,qT2,charge)
    integrand_xib=lambda xib: self.get_M(xia_(xib),xib,x/xia_(xib),z/xib,Q2,muF2,qT2,charge)

    if method=='quad':
      FO = quad(integrand_xia,x+w,1)[0] + quad(integrand_xib,z+w,1)[0]

    elif method=='gauss':
      integrand_xia=np.vectorize(integrand_xia)
      integrand_xib=np.vectorize(integrand_xib)
      FO = fixed_quad(integrand_xia,x+w,1,n=40)[0] + fixed_quad(integrand_xib,z+w,1,n=40)[0]

    if ps=='dxdQ2dzdqT2':
      s=x*y*Q2
      prefactor = D['alphaEM']**2 * self.SC.get_alphaS(muR2) 
      prefactor/= 2*s**2*Q2*x**2
      prefactor*= D['GeV**-2 -> nb'] 
      return prefactor * FO
    else: 
      print 'ps not inplemented'
      return None

class SIDIS_ASY(object):

  def __init__(self,CPDF,CFF):

    self.D={}
    self.setup_params() 
    self.CPDF=CPDF
    self.CFF=CFF
    self.SC=StrongCoupling('one loop')

  def setup_params(self):
    D=self.D
    D['alphaEM']=1/137.035999074
    D['eU2']=(2.0/3.0)**2
    D['eD2']=(-1.0/3.0)**2
    D['GeV**-2 -> nb'] = 0.3894e6
    D['sin2(thetaW)'] = 0.2397
    D['cos2(thetaW)'] = 1-0.2397
    D['CF']=4.0/3.0

  def get_convolution(self,x,splitting,func,method='gauss'):

    D=self.D

    # define integrands

    if splitting=='qq':
      integrand=lambda z: D['CF']/(1-z)*((1+z**2)*func(x/z)-2*func(x))

    elif splitting=='qg':
      integrand=lambda z: 0.5*(1-2*z+2*z*z)*func(x/z)

    elif splitting=='gq':
      integrand=lambda z: D['CF']*(1+(1-z)**2)/z*func(x/z)

    #Z=np.linspace(x,1,100)
    #Y=[integrand(z) for z in Z]
    #py.plot(Z,Y)
    #py.show()
    #sys.exit()

    # compute convolution

    if method=='quad':
      convolution = quad(integrand,x,1)[0]

    elif method=='gauss':
      integrand=np.vectorize(integrand)
      convolution= fixed_quad(integrand,x,1,n=40)[0] 

    # add extra terms from delta

    if splitting=='qq':
      convolution+= D['CF']*(1.5+ 2*np.log(1-x))*func(x)

    return convolution

  def get_PDF_FF(self,x,z,Q2,qT2,muF2,flav,charge,method='gauss'):
    D=self.D

    PDF=lambda x: self.CPDF.get_CPDF(flav,x,muF2)  
    FF =lambda z: self.CFF.get_FF(z,muF2,flav,charge)  
    PDF_g=lambda x: self.CPDF.get_CPDF('g',x,muF2)  
    FF_g =lambda z: self.CFF.get_FF(z,muF2,'g',charge)  

    out = FF(z)*(\
           self.get_convolution(x,'qq',PDF,method=method)\
          +self.get_convolution(x,'qg',PDF_g,method=method))

    out+= PDF(x)*(\
           self.get_convolution(z,'qq',FF,method=method)\
          +self.get_convolution(z,'gq',FF_g,method=method))

    out = 2*D['CF']*FF(z)*PDF(x)*(np.log(Q2/qT2)-1.5)
    return out

  def get_ASY(self,x,y,z,Q2,qT2,muR2,muF2,charge,ps='dxdQ2dzdqT2',method='gauss'):
    D=self.D

    CPDF=self.CPDF
    CFF=self.CFF
    ASY = D['eU2']*self.get_PDF_FF(x,z,Q2,qT2,muF2,'u',charge,method=method)
    ASY+= D['eD2']*self.get_PDF_FF(x,z,Q2,qT2,muF2,'d',charge,method=method)
    ASY+= D['eU2']*self.get_PDF_FF(x,z,Q2,qT2,muF2,'ub',charge,method=method)
    ASY+= D['eD2']*self.get_PDF_FF(x,z,Q2,qT2,muF2,'db',charge,method=method)

    if ps=='dxdQ2dzdqT2':
      s=x*y*Q2
      A1=1+(2/y-1)**2
      prefactor = D['alphaEM']**2 * self.SC.get_alphaS(muR2)*A1
      prefactor/= 2*s**2*x**2*qT2
      prefactor*= D['GeV**-2 -> nb'] 
      return prefactor * ASY

    else: 
      print 'ps not inplemented'
      return None


# testing routines 

def test_SIDIS_L():

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

  sidis=SIDIS_L(CPDF,CFF)
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

def test_SIDIS_FO():

  CPDF=Cteq6PDF('./')
  CFF=FragFuncs('FF/tables/PILO.TAB')

  Q2=15.0
  z=0.5
  x=0.15
  s=300.0
  y=Q2/x/s
  qT=2.0
  charge=1
  muR2=Q2
  muF2=Q2

  sidis=SIDIS_FO(CPDF,CFF)
  #print sidis.get_FO(x,y,z,Q2,qT**2,muR2,muF2,charge,method='quad')
  #print sidis.get_FO(x,y,z,Q2,qT**2,muR2,muF2,charge,method='gauss')
  for i in range(1000):
    print sidis.get_FO(x,y,z,Q2,qT**2,muR2,muF2,charge)

def test_SIDIS_ASY():

  CPDF=Cteq6PDF('./')
  CFF=FragFuncs('FF/tables/PILO.TAB')

  Q2=2.0
  z=0.5
  x=0.005
  s=300.0
  y=Q2/x/s
  charge=0
  muR2=Q2
  muF2=Q2
  qTs=np.linspace(0.1,1.4,40)

  sidis=SIDIS_ASY(CPDF,CFF)
  #print sidis.get_ASY(x,y,z,Q2,qT**2,muR2,muF2,charge,method='quad')
  #print sidis.get_ASY(x,y,z,Q2,qT**2,muR2,muF2,charge,method='gauss')
  #for i in range(1000): print i,sidis.get_ASY(x,y,z,Q2,qT**2,muR2,muF2,charge)

  Y=[sidis.get_ASY(x,y,z,Q2,qT**2,muR2,muF2,charge) for qT in qTs]
  py.plot(qTs,Y)
  py.show()

if __name__== "__main__":

  #test_SIDIS_L()
  #test_SIDIS_FO()
  test_SIDIS_ASY()
 










