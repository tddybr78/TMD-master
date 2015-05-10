#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from cteq6 import cteq6 
    
class Cteq6PDF(object):

  def __init__(self,path2cteq6='./'):
    cteq6.setctq6(1,path2cteq6+'cteq6/')
    self.idx={'g':0,'u':1,'d':2,'s':3,'c':4,'b':5}
    self.idx.update({'ub':-1,'db':-2,'sb':-3,'cb':-4,'bb':-5})

  def get_CPDF(self,flav,x,Q2):
    return cteq6.ctq6pdf(self.idx[flav],x,Q2**0.5)

class ToyPDF(object):

  def get_CPDF(self,flav,x,Q2):

    return x**(-0.1)*(1-x)**3

if __name__== "__main__":

  print dir(cteq6)


