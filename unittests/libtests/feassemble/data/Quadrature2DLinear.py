#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/feassemble/data/Quadrature2DLinear.odb
##
## @brief Python container holding quadrature information for a 2-D
## linear finite-element cell used in testing finite-element C++
## routines.

from pyre.components.Component import Component

import numpy

# ----------------------------------------------------------------------
def N0(p):
  return -0.5*(p[0]+p[1])

def N0p(p):
  return -0.5

def N0q(p):
  return -0.5

def N1(p):
  return 0.5*(1.0+p[0])

def N1p(p):
  return 0.5

def N1q(p):
  return 0.0

def N2(p):
  return 0.5*(1.0+p[1])

def N2p(p):
  return 0.0

def N2q(p):
  return 0.5

# ----------------------------------------------------------------------

# Quadrature2DLinear class
class Quadrature2DLinear(Component):
  """
  Python container holding quadrature information for a 2-D linear
  finite-element cell used in testing finite-element C++ routines.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="quadrature2dlinear"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="quadrature")
    
    self.quadPtsRef = numpy.array( [[-1.0/3.0, -1.0/3.0]], dtype=numpy.float64)
    self.quadWts = numpy.array([2.0], dtype=numpy.float64)
    self.numBasis = 3
    self.numQuadPts = 1
    self.spaceDim = 2
    self.cellDim = 2
    
    return
  

  def calculateBasis(self):
    """
    Calculate basis functions and their derivatives at quadrature points.
    """

    basis = numpy.zeros( (self.numQuadPts, self.numBasis),
                         dtype=numpy.float64)
    basisDeriv = numpy.zeros( (self.numQuadPts, self.numBasis, self.cellDim),
                              dtype=numpy.float64)

    iQuad = 0
    for q in self.quadPtsRef:
      # Basis functions at quadrature points
      basisQ = numpy.array([N0(q), N1(q), N2(q)], dtype=numpy.float64)
      basis[iQuad] = basisQ.reshape( (self.numBasis,) )
      
      # Derivatives of basis functions at quadrature points
      derivQ = numpy.array([[N0p(q), N0q(q)],
                           [N1p(q), N1q(q)],
                           [N2p(q), N2q(q)]], dtype=numpy.float64)      
      basisDeriv[iQuad] = derivQ.reshape((self.numBasis, self.cellDim))

      iQuad += 1
    return (basis, basisDeriv)
    

# FACTORIES ////////////////////////////////////////////////////////////
def quadrature():
  """
  Factory for Quadrature2DLinear.
  """
  return Quadrature2DLinear()


# End of file 
