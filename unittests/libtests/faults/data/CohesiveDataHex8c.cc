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

/* Original mesh
 *
 * Cells are 0-1 and vertices are 2-13.
 *
 *       2,3,4,5 -------- 6,7,8,9 -------- 10,11,12,13
 *
 *                        ^^^^^^^ Vertices forming fault
 *
 * After adding cohesive elements
 *
 * Cells are 0-1,2 and vertices are 3-18.
 *
 *       3,4,5,6 -------- 7,8,9,10 -- 15,16,17,18 -------- 11,12,13,14
 *
 *                        ^^^^^^^^^^^^^^^^^^^^^^ Cohesive element
 *
 */

#include "CohesiveDataHex8c.hh"

const int pylith::faults::CohesiveDataHex8c::_numVertices = 16;

const int pylith::faults::CohesiveDataHex8c::_spaceDim = 3;

const int pylith::faults::CohesiveDataHex8c::_numCells = 3;

const int pylith::faults::CohesiveDataHex8c::_cellDim = 3;

const int pylith::faults::CohesiveDataHex8c::_numCorners[3] = {
  8,
  8,
  8
};

const int pylith::faults::CohesiveDataHex8c::_materialIds[3] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataHex8c::_numGroups = 2;

const int pylith::faults::CohesiveDataHex8c::_groupSizes[2] = {
  8+8+2, 8+8+2 // vertices+edges+faces 
};

const char* pylith::faults::CohesiveDataHex8c::_groupNames[2] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataHex8c::_groupTypes[2] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataHex8c::_filename = 
  "data/hex8c.mesh";

pylith::faults::CohesiveDataHex8c::CohesiveDataHex8c(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numCorners = const_cast<int*>(_numCorners);
  materialIds = const_cast<int*>(_materialIds);
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
  filename = const_cast<char*>(_filename);
} // constructor

pylith::faults::CohesiveDataHex8c::~CohesiveDataHex8c(void)
{}


// End of file
