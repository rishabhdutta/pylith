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

## @file pylith/apps/PyLithApp.py
##
## @brief Python PyLith application

from PetscApplication import PetscApplication

# PyLithApp class
class PyLithApp(PetscApplication):
  """
  Python PyLithApp application.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscApplication.Inventory):
    """
    Python object for managing PyLithApp facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing PyLithApp facilities and properties.
    ##
    ## \b Properties
    ## @li \b initialize_only Stop simulation after initializing problem.
    ##
    ## \b Facilities
    ## @li \b mesher Generates or imports the computational mesh.
    ## @li \b problem Computational problem to solve
    ## @li \b petsc Manager for PETSc options

    import pyre.inventory

    initializeOnly = pyre.inventory.bool("initialize_only", default=False)
    initializeOnly.meta['tip'] = "Stop simulation after initializing problem."

    from pylith.topology.MeshImporter import MeshImporter
    mesher = pyre.inventory.facility("mesh_generator", family="mesh_generator",
                                     factory=MeshImporter)
    mesher.meta['tip'] = "Generates or imports the computational mesh."

    from pylith.problems.TimeDependent import TimeDependent
    problem = pyre.inventory.facility("problem", family="problem",
                                      factory=TimeDependent)
    problem.meta['tip'] = "Computational problem to solve."

    from pylith.perf.MemoryLogger import MemoryLogger
    perfLogger = pyre.inventory.facility("perf_logger", family="perf_logger",
                                         factory=MemoryLogger)
    perfLogger.meta['tip'] = "Performance and memory logging."

    
    typos = pyre.inventory.str("typos", default="pedantic",
                               validator=pyre.inventory.choice(['relaxed', 'strict', 'pedantic']))
    typos.meta['tip'] = "Specifies the handling of unknown properties and " \
        "facilities"
    
    pdbOn = pyre.inventory.bool("start_python_debugger", default=False)
    pdbOn.meta['tip'] = "Start python debugger at beginning of main()."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="pylithapp"):
    """
    Constructor.
    """
    PetscApplication.__init__(self, name)
    self._loggingPrefix = "PyLith "
    return


  def main(self, *args, **kwds):
    """
    Run the application.
    """
    if self.pdbOn:
          import pdb
          pdb.set_trace()
        
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()
    if 0 == comm.rank:
      self._info.log("Running on %d process(es)." % comm.size)

    from pylith.utils.profiling import resourceUsageString    
    self._debug.log(resourceUsageString())

    self._setupLogging()

    # Create mesh (adjust to account for interfaces (faults) if necessary)
    self._eventLogger.stagePush("Meshing")
    interfaces = None
    if "interfaces" in dir(self.problem):
      interfaces = self.problem.interfaces.components()
    mesh = self.mesher.create(self.problem.normalizer, interfaces)
    del interfaces
    del self.mesher
    self._debug.log(resourceUsageString())
    self._eventLogger.stagePop()

    # Setup problem, verify configuration, and then initialize
    self._eventLogger.stagePush("Setup")
    self.problem.preinitialize(mesh)
    self._debug.log(resourceUsageString())

    self.problem.verifyConfiguration()

    self.problem.initialize()
    self._debug.log(resourceUsageString())

    self._eventLogger.stagePop()

    # If initializing only, stop before running problem
    if self.initializeOnly:
      return

    # Run problem
    self.problem.run(self)
    self._debug.log(resourceUsageString())

    # Cleanup
    self._eventLogger.stagePush("Finalize")
    self.problem.finalize()
    self._eventLogger.stagePop()

    self.perfLogger.logMesh('Mesh', mesh)
    self.compilePerformanceLog()
    if self.perfLogger.verbose:
      self.perfLogger.show()

    return


  def version(self):
    def getPyPkgVer(name):
      import os
      m = None
      location = None
      version = None
      try:
        m = __import__(name)
        location = os.path.split(m.__file__)[0]
        version = m.__version__
      except ImportError:
        version = "not found"
        location = "--"
      except AttributeError:
        if version is None:
          version = "unknown"
        if location is None:
          location = "unknown"
      return (version, location)
    
    import sys
    import platform
    pythonVersion = platform.python_version()
    pythonCompiler = platform.python_compiler()
    uname = platform.uname()
    unameStr = " ".join((uname[0], uname[2], uname[4]))
    msg = "Running PyLith on %s.\n" % unameStr

    import pylith.utils.utils as utils
    # PyLith version information
    v = utils.PylithVersion()
    if v.isRelease():
        msg += "    Release v%s.\n" % (v.version(),)
    else:
        msg += "    Configured on %s, GIT branch: %s, revision: %s, hash: %s.\n" % (v.gitDate(), v.gitBranch(), v.gitRevision(), v.gitHash(),)
    msg += "\n"
        
    # PETSc
    v = utils.PetscVersion()
    if v.isRelease():
        msg += "    PETSc release v%s.\n" % (v.version(),)
    else:
        msg += "    PETSc configured on %s, GIT branch: %s, revision: %s.\n" % (v.gitDate(), v.gitBranch(), v.gitRevision(),)
    msg += "        PETSC_DIR: %s, PETSC_ARCH: %s\n" % (v.petscDir(), v.petscArch(),)
    msg += "\n"

    # Other dependencies
    v = utils.DependenciesVersion()
    msg += "    MPI standard: %s, implementation: %s, version: %s.\n" % (v.mpiStandard(), v.mpiImplementation(), v.mpiVersion())
    msg += "    HDF5 version: %s.\n" % (v.hdf5Version())
    msg += "    NetCDF4 version: %s.\n" % (v.netcdfVersion())
    msg += "\n"
    
    # Spatialdata
    import spatialdata.utils.utils as utils
    v = utils.SpatialdataVersion()
    if v.isRelease():
        msg += "    Spatialdata release v%s.\n" % (v.version(),)
    else:
        msg += "    Spatialdata configured on %s, GIT branch: %s, revision: %s.\n" % (v.gitDate(), v.gitBranch(), v.gitRevision(),)
    msg += "    Proj.4 version: %s.\n" % (v.projVersion(),)
    msg += "\n"

    # Python
    msg += "    Python %s compiled with %s from %s.\n" % (pythonVersion, pythonCompiler, sys.executable)

    pkgs = ("numpy","spatialdata","FIAT","h5py","netCDF4","pyre")
    for pkg in pkgs:
      ver,loc = getPyPkgVer(pkg)
      msg += "        %s %s from %s.\n" % (pkg, ver, loc)
    print(msg)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    PetscApplication._configure(self)
    self.initializeOnly = self.inventory.initializeOnly
    self.mesher = self.inventory.mesher
    self.problem = self.inventory.problem
    self.perfLogger = self.inventory.perfLogger
    self.typos = self.inventory.typos
    self.pdbOn = self.inventory.pdbOn

    import journal
    self._debug = journal.debug(self.name)
    return

  def _setupLogging(self):
    """
    Setup event logging.
    """
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("PyLith")
    logger.initialize()

    self._eventLogger = logger
    return
  

# ======================================================================
# Local version of InfoApp that only configures itself. Workaround for
# not adding --help-all (or similar property) in Pyre Application.
#
# Note: We want --help-all to display settings before launching a
# parallel job.
class InfoApp(PyLithApp):
  def __init__(self, args, name="pylithapp"):
    """
    Constructor.
    """
    PyLithApp.__init__(self, name)
    self.pylithargs = args
    return

  def onLoginNode(self, *args, **kwds):
    """
    Instead of scheduling job, do nothing.
    """
    return

  def getArgv(self, *args, **kwds):
    """
    Prevent PyLith from getting all of the command line arguments. Use
    only the ones relevant to PyLith which are specified in the arg to
    the constructor.
    """
    argv = kwds.get('argv')
    if argv is None:
      argv = self.pylithargs
    else:
      self.arg0 = argv[0]
      self._requires = kwds.get('requires')
      argv = argv[1:]
    return argv
  



# End of file 
