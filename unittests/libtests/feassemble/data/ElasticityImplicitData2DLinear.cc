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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application elasticityapp.

#include "ElasticityImplicitData2DLinear.hh"

const int pylith::feassemble::ElasticityImplicitData2DLinear::_spaceDim = 2;

const int pylith::feassemble::ElasticityImplicitData2DLinear::_cellDim = 2;

const int pylith::feassemble::ElasticityImplicitData2DLinear::_numVertices = 3;

const int pylith::feassemble::ElasticityImplicitData2DLinear::_numCells = 1;

const int pylith::feassemble::ElasticityImplicitData2DLinear::_numBasis = 3;

const int pylith::feassemble::ElasticityImplicitData2DLinear::_numQuadPts = 1;

const char* pylith::feassemble::ElasticityImplicitData2DLinear::_matType = "ElasticPlaneStrain";

const char* pylith::feassemble::ElasticityImplicitData2DLinear::_matDBFilename = "data/elasticplanestrain.spatialdb";

const int pylith::feassemble::ElasticityImplicitData2DLinear::_matId = 0;

const char* pylith::feassemble::ElasticityImplicitData2DLinear::_matLabel = "elastic strain 2-D";

const double pylith::feassemble::ElasticityImplicitData2DLinear::_dt =   1.00000000e-02;

const double pylith::feassemble::ElasticityImplicitData2DLinear::_gravityVec[] = {
  0.00000000e+00, -1.00000000e+08,
};

const double pylith::feassemble::ElasticityImplicitData2DLinear::_vertices[] = {
  2.00000000e-01, -4.00000000e-01,
  3.00000000e-01,  5.00000000e-01,
 -1.00000000e+00, -2.00000000e-01,
};

const int pylith::feassemble::ElasticityImplicitData2DLinear::_cells[] = {
0,1,2,
};

const double pylith::feassemble::ElasticityImplicitData2DLinear::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::ElasticityImplicitData2DLinear::_quadPts[] = {
 -3.33333333e-01, -3.33333333e-01,
};

const double pylith::feassemble::ElasticityImplicitData2DLinear::_quadWts[] = {
  2.00000000e+00,
};

const double pylith::feassemble::ElasticityImplicitData2DLinear::_basis[] = {
  3.33333333e-01,  3.33333333e-01,
  3.33333333e-01,};

const double pylith::feassemble::ElasticityImplicitData2DLinear::_basisDerivRef[] = {
 -5.00000000e-01, -5.00000000e-01,
  5.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  5.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitData2DLinear::_fieldTIncr[] = {
  1.30000000e+00, -9.00000000e-01,
  1.40000000e+00,  1.50000000e+00,
  5.00000000e-01, -9.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitData2DLinear::_fieldT[] = {
  1.60000000e+00, -8.00000000e-01,
  9.00000000e-01,  7.00000000e-01,
 -2.00000000e-01, -1.10000000e+00,
};

const double pylith::feassemble::ElasticityImplicitData2DLinear::_fieldTmdt[] = {
  8.00000000e-01,  1.00000000e-01,
  5.00000000e-01,  3.00000000e-01,
 -1.00000000e-01, -6.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitData2DLinear::_valsResidual[] = {
 -1.04842898e+11,  2.84328125e+11,
 -3.14863636e+10, -2.63281250e+11,
  1.36329261e+11, -2.10468750e+10,
};

const double pylith::feassemble::ElasticityImplicitData2DLinear::_valsJacobian[] = {
  4.35710227e+10, -2.45596591e+10,
 -1.59886364e+10,  7.35795455e+09,
 -2.75823864e+10,  1.72017045e+10,
 -2.45596591e+10,  7.59573864e+10,
  8.29545455e+09, -6.18693182e+10,
  1.62642045e+10, -1.40880682e+10,
 -1.59886364e+10,  8.29545455e+09,
  2.16818182e+10,  6.47727273e+09,
 -5.69318182e+09, -1.47727273e+10,
  7.35795455e+09, -6.18693182e+10,
  6.47727273e+09,  5.94659091e+10,
 -1.38352273e+10,  2.40340909e+09,
 -2.75823864e+10,  1.62642045e+10,
 -5.69318182e+09, -1.38352273e+10,
  3.32755682e+10, -2.42897727e+09,
  1.72017045e+10, -1.40880682e+10,
 -1.47727273e+10,  2.40340909e+09,
 -2.42897727e+09,  1.16846591e+10,
};

pylith::feassemble::ElasticityImplicitData2DLinear::ElasticityImplicitData2DLinear(void)
{ // constructor
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numVertices = _numVertices;
  numCells = _numCells;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  matType = const_cast<char*>(_matType);
  matDBFilename = const_cast<char*>(_matDBFilename);
  matId = _matId;
  matLabel = const_cast<char*>(_matLabel);
  dt = _dt;
  gravityVec = const_cast<double*>(_gravityVec);
  vertices = const_cast<double*>(_vertices);
  cells = const_cast<int*>(_cells);
  verticesRef = const_cast<double*>(_verticesRef);
  quadPts = const_cast<double*>(_quadPts);
  quadWts = const_cast<double*>(_quadWts);
  basis = const_cast<double*>(_basis);
  basisDerivRef = const_cast<double*>(_basisDerivRef);
  fieldTIncr = const_cast<double*>(_fieldTIncr);
  fieldT = const_cast<double*>(_fieldT);
  fieldTmdt = const_cast<double*>(_fieldTmdt);
  valsResidual = const_cast<double*>(_valsResidual);
  valsJacobian = const_cast<double*>(_valsJacobian);
} // constructor

pylith::feassemble::ElasticityImplicitData2DLinear::~ElasticityImplicitData2DLinear(void)
{}


// End of file
