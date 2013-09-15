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

#include <portinfo>

#include "FieldOps.hh" // implementation of class methods

#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR

// ----------------------------------------------------------------------
void
pylith::topology::FieldOps::createFiniteElement(std::string name, PetscInt dim, PetscInt numComponents, PetscInt order, PetscFE *fe)
{ // createFiniteElement

  PetscFE         fem;
  PetscQuadrature q;
  DM              K;
  PetscSpace      P;
  PetscDualSpace  Q;
  std::string     prefix = name + "_";
  PetscErrorCode  err;

  PYLITH_METHOD_BEGIN;
  /* Create space */
  err = PetscSpaceCreate(PETSC_COMM_SELF, &P);PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetOptionsPrefix((PetscObject) P, prefix.c_str());PYLITH_CHECK_ERROR(err);
  err = PetscSpaceSetFromOptions(P);PYLITH_CHECK_ERROR(err);
  err = PetscSpacePolynomialSetNumVariables(P, dim);PYLITH_CHECK_ERROR(err);
  err = PetscSpaceSetOrder(P, order);PYLITH_CHECK_ERROR(err);
  err = PetscSpaceSetUp(P);PYLITH_CHECK_ERROR(err);
  /* Create dual space */
  err = PetscDualSpaceCreate(PETSC_COMM_SELF, &Q);PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetOptionsPrefix((PetscObject) Q, prefix.c_str());PYLITH_CHECK_ERROR(err);
  err = PetscDualSpaceCreateReferenceCell(Q, dim, PETSC_TRUE, &K);PYLITH_CHECK_ERROR(err);
  err = PetscDualSpaceSetDM(Q, K);PYLITH_CHECK_ERROR(err);
  err = DMDestroy(&K);PYLITH_CHECK_ERROR(err);
  err = PetscDualSpaceSetOrder(Q, order);PYLITH_CHECK_ERROR(err);
  err = PetscDualSpaceSetFromOptions(Q);PYLITH_CHECK_ERROR(err);
  err = PetscDualSpaceSetUp(Q);PYLITH_CHECK_ERROR(err);
  /* Create element */
  err = PetscFECreate(PETSC_COMM_SELF, &fem);PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetOptionsPrefix((PetscObject) fem, prefix.c_str());PYLITH_CHECK_ERROR(err);
  err = PetscFESetFromOptions(fem);PYLITH_CHECK_ERROR(err);
  err = PetscFESetBasisSpace(fem, P);PYLITH_CHECK_ERROR(err);
  err = PetscFESetDualSpace(fem, Q);PYLITH_CHECK_ERROR(err);
  err = PetscFESetNumComponents(fem, numComponents);PYLITH_CHECK_ERROR(err);
  err = PetscSpaceDestroy(&P);PYLITH_CHECK_ERROR(err);
  err = PetscDualSpaceDestroy(&Q);PYLITH_CHECK_ERROR(err);
  /* Create quadrature */
  err = PetscDTGaussJacobiQuadrature(dim, order, -1.0, 1.0, &q);PYLITH_CHECK_ERROR(err);
  err = PetscFESetQuadrature(fem, q);PYLITH_CHECK_ERROR(err);
  *fe = fem;
  PYLITH_METHOD_END;
} // createFiniteElement

// End of file 
