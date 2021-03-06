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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#if !defined(pylith_meshio_geomdatatet3d_hh)
#define pylith_meshio_geomdatatet3d_hh

#include "CellGeomData.hh"

namespace pylith {
  namespace feassemble {
     class GeomDataTet3D;
  } // feassemble
} // pylith

class pylith::feassemble::GeomDataTet3D : public CellGeomData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  GeomDataTet3D(void);

  /// Destructor
  ~GeomDataTet3D(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const int _cellDim; ///< Number of dimensions associated with cell
  static const int _spaceDim; ///< Number of dimensions in vertex coordinates
  static const int _numCorners; ///< Number of vertices in cell

  static const int _numLocs; ///< Number of locations for computing Jacobian

  static const PylithScalar _gravityVec[]; ///< Constant gravity vector
  static const PylithScalar _vertices[]; ///< Coordinates of cell's vertices
  static const PylithScalar _locations[]; ///< Locations to compute Jacobian
  static const PylithScalar _jacobian[]; ///< Jacobian at locations
  static const PylithScalar _jacobianDet[]; ///< Determinant of Jacobian at locations

};

#endif // pylith_meshio_geomdatatet3d_hh

// End of file
