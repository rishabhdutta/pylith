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
# Copyright (c) 2010-2013 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/feassemble/Discretization.py
##
## @brief Python object for managing basis functions and quadrature
## parameters.
##
## Factory: discretization.

from pylith.utils.PetscComponent import PetscComponent

# Discretization class
class Discretization(PetscComponent):
  """
  Python object for managing basis functions and quadrature parameters.

  Factory: discretization.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """Python object for managing Discretization facilities and properties."""

    ## @class Inventory
    ## Python object for managing Discretization facilities and properties.
    ##
    ## \b Properties
    ## @li \b basis_order Order of basis functions.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    basisOrder = pyre.inventory.int("basis_order", default=1, validator=pyre.inventory.greaterEqual(0))
    basisOrder.meta['tip'] = "Order of basis functions."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="discretization"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    try:
      PetscComponent._configure(self)
      self.basisOrder = self.inventory.basisOrder

    except ValueError, err:
      aliases = ", ".join(self.aliases)
      raise ValueError("Error while configuring discretization "
                       "(%s):\n%s" % (aliases, err.message))
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def discretization():
  """
  Factory associated with Discretization.
  """
  return Discretization()


# End of file 
