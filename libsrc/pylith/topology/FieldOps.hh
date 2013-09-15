// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/topology/FieldOps.hh
 *
 * @brief Utility functions for fields.
 */

#if !defined(pylith_topology_fieldops_hh)
#define pylith_topology_fieldops_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations
#include "pylith/utils/petscfwd.h" // HASA PetscVec

#include <string> // USES std::string

// FieldOps ----------------------------------------------------------------
/** @brief Utility functions for fields.
 *
 * Utility functions for fields
 */
class pylith::topology::FieldOps
{ // FieldOps
public:
  static void createFiniteElement(std::string name, PetscInt dim, PetscInt numComponents, PetscInt order, PetscFE *fe);
}; // FieldOps

#endif // pylith_topology_fieldops_hh

// End of file 
