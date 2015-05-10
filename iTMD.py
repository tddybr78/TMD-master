#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from IPython.html.widgets import interact, interactive, fixed
from IPython.html import widgets
from IPython.display import clear_output, display, HTML
from CollinearPDF import Cteq6PDF,FakePDF

class InteractiveTMD(object):

  def __init__(self,root):
    self.root=root

  def select_CPDFS(self):

    def f(**kwargs):
      if kwargs['pdfset']=='cteq6': self.CPDF=Cteq6PDF(self.root)
      if kwargs['pdfset']=='fake': self.CPDF=FakePDF()

    interact(f,pdfset={'cteq6':'cteq6','fake':'fake'})
 
    


