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

## @file pylith/meshio/ElasticityFields.py
##
## @brief Python container with fields in solution for elasticty problems.

from pylith.utils.PetscComponent import PetscComponent

# ElasticityFields class
class ElasticityFields(PetscComponent):
  """
  Python container with fields in solution for elasticity problems.

  Factory: object_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing ElasticityFields facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing ElasticityFields facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b displacement Displacement field discretization.

    import pyre.inventory

    from Discretization import Discretization
    displacement = pyre.inventory.facility("displacement", family="discretization", factory=Discretization)
    displacement.meta['tip'] = "Displacement field discretization."


  def initialze(self, dimension):
      """
      Initialize fields.
      """
      displacement.dimension = dimension
      return
  
  
  def components(self):
      """
      Get fields in desired order.
      """
      return [self.displacement]


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticityfields"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="discretization")
    return


# End of file 
