#!/usr/bin/env python
# Author: Nobuo Sato (nsato@jlab.org)
import sys

# This is a simple example that shows how to 
# use the TMD python modules

# we need to tell python kwhere the TMD modules
# located
sys.path.insert(1,'../../')

# now lets load the TMD module
from TMDPDF import TMDPDF

# lets create an instance (object) of the TMDPDF class.
# we will call it TMD

TMD = TMDPDF()

# a function in the class called "get_PDF" is 
# is a routine that computes TMD PDF in qT space.
# we can acess sinng the object that we have created 
# (e.g: TMD.name_of_func(...))

# first lets specify some parameters
x=0.5
mu2=5
xiF=mu2
qT=2

print TMD.get_PDF(x,mu2,xiF,qT)

# end of example :)
