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

## @file tests/2d/quad4/friction_compression_soln.py
##
## @brief Analytical solution to compression problem with friction.

import numpy

# Physical properties
p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density*p_vs**2
p_lambda = p_density*p_vp**2 - 2*p_mu

# Uniform stress field (plane strain)
sxx = 0.0
sxy = 1.0e+6
syy = 0.0
szz = p_lambda/(2*p_lambda+2*p_mu)*(sxx+syy)

# Uniform strain field
exx = 1.0/(2*p_mu) * (sxx - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz))
eyy = 1.0/(2*p_mu) * (syy - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz))
ezz = 1.0/(2*p_mu) * (szz - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz))

exy = 1.0/(2*p_mu) * (sxy)

#print exx,eyy,exy,ezz,szz
#print -exx*p_lambda/(p_lambda+2*p_mu)

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
  """
  Analytical solution to friction_compression problem.
  """

  def __init__(self):
    return


  def displacement(self, locs, nlocsO):
    """
    Compute displacement field at locations.
    """
    (nlocs, dim) = locs.shape

    disp = numpy.zeros( (1, nlocs, 2), dtype=numpy.float64)
    disp[0,:,1] = 2*exy*(locs[:,0]+max(abs(locs[:,0])))
    return disp


  def strain(self, locs):
    """
    Compute strain field at locations.
    """
    (npts, dim) = locs.shape
    strain = numpy.zeros( (1, npts, 3), dtype=numpy.float64)
    strain[0,:,0] = exx
    strain[0,:,1] = eyy
    strain[0,:,2] = exy
    return strain
  

  def stress(self, locs):
    """
    Compute stress field at locations.
    """
    (npts, dim) = locs.shape
    stress = numpy.zeros( (1, npts, 3), dtype=numpy.float64)
    stress[0,:,0] = sxx
    stress[0,:,1] = syy
    stress[0,:,2] = sxy
    return stress


# End of file 
