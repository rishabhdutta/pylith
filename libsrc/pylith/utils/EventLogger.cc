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

#include <portinfo>

#include "EventLogger.hh" // Implementation of class methods

#include "error.h" // USES PYLITH_METHOD_BEGIN/END

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::utils::EventLogger::EventLogger(void) :
  _className(""),
  _classId(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::utils::EventLogger::~EventLogger(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Setup logging class.
void
pylith::utils::EventLogger::initialize(void)
{ // initialize
  PYLITH_METHOD_BEGIN;

  if (_className == "")
    throw std::logic_error("Must set logging class name before "
			   "initializaing EventLogger.");
  
  _events.clear();
  PetscErrorCode err = PetscClassIdRegister(_className.c_str(), &_classId);
  if (err) {
    std::ostringstream msg;
    msg << "Could not register logging class '" << _className << "'.";
    throw std::runtime_error(msg.str());
  } // if
  assert(_classId);

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Register event.
int
pylith::utils::EventLogger::registerEvent(const char* name)
{ // registerEvent
  PYLITH_METHOD_BEGIN;

  assert(_classId);
  int id = 0;
  PetscErrorCode err = PetscLogEventRegister(name, _classId, &id);
  if (err) {
    std::ostringstream msg;
    msg << "Could not register logging event '" << name
	<< "' for logging class '" << _className << "'.";
    throw std::runtime_error(msg.str());
  } // if  
  _events[name] = id;
  PYLITH_METHOD_RETURN(id);
} // registerEvent

// ----------------------------------------------------------------------
// Get event identifier.
int
pylith::utils::EventLogger::eventId(const char* name)
{ // eventId
  PYLITH_METHOD_BEGIN;

  map_event_type::iterator iter = _events.find(name);
  if (iter == _events.end()) {
    std::ostringstream msg;
    msg << "Could not find logging event '" << name
	<< "' in logging class '" << _className << "'.";
    throw std::runtime_error(msg.str());
  } // if

  PYLITH_METHOD_RETURN(iter->second);
} // eventId

// ----------------------------------------------------------------------
// Register stage.
int
pylith::utils::EventLogger::registerStage(const char* name)
{ // registerStage
  PYLITH_METHOD_BEGIN;

  assert(_classId);
  int id = 0;
  PetscErrorCode err = PetscLogStageRegister(name, &id);
  if (err) {
    std::ostringstream msg;
    msg << "Could not register logging stage '" << name << "'.";
    throw std::runtime_error(msg.str());
  } // if  
  _stages[name] = id;

  PYLITH_METHOD_RETURN(id);
} // registerStage

// ----------------------------------------------------------------------
// Get stage identifier.
int
pylith::utils::EventLogger::stageId(const char* name)
{ // stageId
  PYLITH_METHOD_BEGIN;

  map_event_type::iterator iter = _stages.find(name);
  if (iter == _stages.end()) {
    registerStage(name);
    iter = _stages.find(name);
    assert(iter != _stages.end());
  } // if

  PYLITH_METHOD_RETURN(iter->second);
} // stagesId


// End of file 
