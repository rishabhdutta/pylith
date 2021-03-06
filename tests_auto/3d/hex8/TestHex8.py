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

## @file tests/3d/hex8/TestHex8.py
##
## @brief Generic tests for problems using 3-D mesh.

import unittest
import numpy

from pylith.tests import has_h5py

class TestHex8(unittest.TestCase):
  """
  Generic tests for problems using 3-D hex8 mesh.
  """

  def setUp(self):
    """
    Setup for tests.
    """
    self.mesh = {'ncells-elastic': 180,
                 'ncells-viscoelastic': 180,
                 'ncorners': 8,
                 'nvertices': 550,
                 'spaceDim': 3,
                 'tensorSize': 6}

    if has_h5py():
      self.checkResults = True
    else:
      self.checkResults = False
    return


  def test_elastic_info(self):
    """
    Check elastic info.
    """
    if not self.checkResults:
      return

    for material in ["elastic", "viscoelastic"]:
      ncells= self.mesh['ncells-%s' % material]
      self.mesh['ncells'] = ncells

      filename = "%s-%s_info.h5" % (self.outputRoot, material)
      from axialdisp_soln import p_mu,p_lambda,p_density

      propMu =  p_mu*numpy.ones( (1, ncells, 1), dtype=numpy.float64)
      propLambda = p_lambda*numpy.ones( (1, ncells, 1), dtype=numpy.float64)
      propDensity = p_density*numpy.ones( (1, ncells, 1), dtype=numpy.float64)

      properties = {'mu': propMu,
                    'lambda': propLambda,
                    'density': propDensity}

      from pylith.tests.PhysicalProperties import check_properties
      check_properties(self, filename, self.mesh, properties)

    return


  def test_soln(self):
    """
    Check solution (displacement) field.
    """
    if not self.checkResults:
      return

    filename = "%s.h5" % self.outputRoot
    from pylith.tests.Solution import check_displacements
    check_displacements(self, filename, self.mesh)

    return


  def test_elastic_statevars(self):
    """
    Check elastic state variables.
    """
    if not self.checkResults:
      return


    for material in ["elastic", "viscoelastic"]:
      ncells= self.mesh['ncells-%s' % material]
      self.mesh['ncells'] = ncells

      filename = "%s-%s.h5" % (self.outputRoot, material)

      from pylith.tests.StateVariables import check_state_variables
      stateVars = ["total_strain", "stress"]
      check_state_variables(self, filename, self.mesh, stateVars)

    return


# End of file
