#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from scipy.interpolate import interp2d
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.integrate import quad,dblquad,fixed_quad
#from line_profiler import LineProfiler

class FragFuncs(object):

  def __init__(self,fname):

    self.D={}
    self.load_table(fname)

  def load_table(self,fname):

    # read file
    F=open(fname,'r')
    L=F.readlines()
    F.close()

    # convert loaded file to floats
    L=[l.strip() for l in L]
    L=[[float(x) for x in l.split()] for l in L]

    # extract parts
    Q2=np.copy(L[0])
    X=np.copy(L[1])
    L=np.array(L[2:])

    # get number of partons in table
    npartons=len(L[0])

    # empy array for FF values
    FF=np.zeros((npartons,Q2.size,X.size))

    # fill array
    cnt=0
    for iX in range(X.size-1):
      for iQ2 in range(Q2.size):
        for iparton in range(npartons):
          if any([iparton==k for k in [0,1,2,6,7,8]]) : 
            factor=(1-X[iX])**4 * X[iX]**0.5
          elif any([iparton==k for k in [3,4]]) : 
            factor=(1-X[iX])**7 * X[iX]**0.3
          elif iparton==5: 
            factor=(1-X[iX])**4 * X[iX]**0.3
          FF[iparton,iQ2,iX]=L[cnt,iparton]/factor
        cnt+=1

    LX=np.log(X)
    LQ2=np.log(Q2)

    D = self.D
    D['UTOT']=RectBivariateSpline(LQ2,LX,FF[0],kx=1,ky=1)
    D['DTOT']=RectBivariateSpline(LQ2,LX,FF[1],kx=1,ky=1)
    D['STOT']=RectBivariateSpline(LQ2,LX,FF[2],kx=1,ky=1)
    D['CTOT']=RectBivariateSpline(LQ2,LX,FF[3],kx=1,ky=1)
    D['BTOT']=RectBivariateSpline(LQ2,LX,FF[4],kx=1,ky=1)
    D['G']   =RectBivariateSpline(LQ2,LX,FF[5],kx=1,ky=1)
    D['UVAL']=RectBivariateSpline(LQ2,LX,FF[6],kx=1,ky=1)
    D['DVAL']=RectBivariateSpline(LQ2,LX,FF[7],kx=1,ky=1)
    D['SVAL']=RectBivariateSpline(LQ2,LX,FF[8],kx=1,ky=1)

    #D['Up'] =RectBivariateSpline(LQ2,LX,0.5*(FF[0]+FF[6]))
    #D['UBp']=RectBivariateSpline(LQ2,LX,0.5*(FF[0]-FF[6]))
    #D['Dp'] =RectBivariateSpline(LQ2,LX,0.5*(FF[1]+FF[7]))
    #D['DBp']=RectBivariateSpline(LQ2,LX,0.5*(FF[1]-FF[7]))
    #D['Sp'] =RectBivariateSpline(LQ2,LX,0.5*(FF[2]+FF[8]))
    #D['SBp']=RectBivariateSpline(LQ2,LX,0.5*(FF[2]-FF[8]))
    #D['Cp'] =RectBivariateSpline(LQ2,LX,0.5*FF[3])
    #D['Bp'] =RectBivariateSpline(LQ2,LX,0.5*FF[4])
    #D['G']  =RectBivariateSpline(LQ2,LX,FF[5])

    D['Up'] = lambda lQ2,lx: 0.5*(D['UTOT'](lQ2,lx)+D['UVAL'](lQ2,lx))
    D['UBp']= lambda lQ2,lx: 0.5*(D['UTOT'](lQ2,lx)-D['UVAL'](lQ2,lx))
    D['Dp'] = lambda lQ2,lx: 0.5*(D['DTOT'](lQ2,lx)+D['DVAL'](lQ2,lx))
    D['DBp']= lambda lQ2,lx: 0.5*(D['DTOT'](lQ2,lx)-D['DVAL'](lQ2,lx))
    D['Sp'] = lambda lQ2,lx: 0.5*(D['STOT'](lQ2,lx)+D['SVAL'](lQ2,lx))
    D['SBp']= lambda lQ2,lx: 0.5*(D['STOT'](lQ2,lx)-D['SVAL'](lQ2,lx))
    D['Cp'] = lambda lQ2,lx: 0.5*D['CTOT'](lQ2,lx)
    D['Bp'] = lambda lQ2,lx: 0.5*D['BTOT'](lQ2,lx)

  def get_xg(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    return self.D['G'](lQ2,lx)[0,0]*(1-x)**4*x**0.3

  def get_xu(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    if charge==1:
      return D['Up'](lQ2,lx)[0,0]*(1-x)**4*x**0.5
    elif charge==-1:
      return D['UBp'](lQ2,lx)[0,0]*(1-x)**4*x**0.5
    elif charge==0:
      return 0.5 *(D['Up'](lQ2,lx)[0,0]+D['UBp'](lQ2,lx)[0,0])*(1-x)**4*x**0.5

  def get_xub(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    if charge==1:
      return D['UBp'](lQ2,lx)[0,0]*(1-x)**4*x**0.5
    elif charge==-1:
      return D['Up'](lQ2,lx)[0,0]*(1-x)**4*x**0.5
    elif charge==0:
      return 0.5*(D['Up'](lQ2,lx)[0,0]+D['UBp'](lQ2,lx)[0,0])*(1-x)**4*x**0.5

  def get_xd(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    if charge==1:
      return D['Dp'](lQ2,lx)[0,0]*(1-x)**4*x**0.5
    elif charge==-1:
      return D['DBp'](lQ2,lx)[0,0]*(1-x)**4*x**0.5
    elif charge==0:
      return 0.5*(D['Dp'](lQ2,lx)[0,0]+D['DBp'](lQ2,lx)[0,0])*(1-x)**4*x**0.5

  def get_xdb(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    if charge==1:
      return D['DBp'](lQ2,lx)[0,0]*(1-x)**4*x**0.5
    elif charge==-1:
      return D['Dp'](lQ2,lx)[0,0]*(1-x)**4*x**0.5
    elif charge==0:
      return 0.5*(D['Dp'](lQ2,lx)[0,0]+D['DBp'](lQ2,lx)[0,0])*(1-x)**4*x**0.5

  def get_xs(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    if charge==1:
      return D['Sp'](lQ2,lx)[0,0]*(1-x)**4*x**0.5
    elif charge==-1:
      return D['SBp'](lQ2,lx)[0,0]*(1-x)**4*x**0.5
    elif charge==0:
      return 0.5*(D['Sp'](lQ2,lx)[0,0]+D['SBp'](lQ2,lx)[0,0])*(1-x)**4*x**0.5

  def get_xsb(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    if charge==1:
      return D['SBp'](lQ2,lx)[0,0]*(1-x)**4*x**0.5
    elif charge==-1:
      return D['Sp'](lQ2,lx)[0,0]*(1-x)**4*x**0.5
    elif charge==0:
      return 0.5*(D['Sp'](lQ2,lx)[0,0]+D['SBp'](lQ2,lx)[0,0])*(1-x)**4*x**0.5

  def get_xc(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    return D['Cp'](lQ2,lx)[0,0]*(1-x)**7*x**0.3

  def get_xcb(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    return D['Cp'](lQ2,lx)[0,0]*(1-x)**7*x**0.3

  def get_xb(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    return D['Bp'](lQ2,lx)[0,0]*(1-x)**7*x**0.3

  def get_xbb(self,x,Q2,charge):
    lx=np.log(x)
    lQ2=np.log(Q2)
    D=self.D
    return D['Bp'](lQ2,lx)[0,0]*(1-x)**7*x**0.3

  def get_FF(self,x,Q2,flav,charge):
    if   flav=='g':  return self.get_xg(x,Q2,charge)/x
    elif flav=='u':  return self.get_xu(x,Q2,charge)/x
    elif flav=='ub': return self.get_xub(x,Q2,charge)/x
    elif flav=='d':  return self.get_xd(x,Q2,charge)/x
    elif flav=='db': return self.get_xdb(x,Q2,charge)/x
    elif flav=='s':  return self.get_xs(x,Q2,charge)/x
    elif flav=='sb': return self.get_xsb(x,Q2,charge)/x
    elif flav=='c':  return self.get_xc(x,Q2,charge)/x
    elif flav=='cb': return self.get_xcb(x,Q2,charge)/x
    elif flav=='b':  return self.get_xb(x,Q2,charge)/x
    elif flav=='bb': return self.get_xbb(x,Q2,charge)/x


if __name__== "__main__":

  import fDSS

  PI=FragFuncs('tables/PILO.TAB')
  ic=0
  print PI.get_xg(0.6,10.0,ic)
  print fDSS.fdss(1,ic,0,0.6,10.0)[8]

























