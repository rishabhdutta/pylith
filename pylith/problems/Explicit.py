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

## @file pylith/problems/Explicit.py
##
## @brief Python Explicit object for solving equations using an
## explicit formulation with a lumped Jacobian that is stored as a
## Field.
##
## Factory: pde_formulation

from Formulation import Formulation
from problems import Explicit as ModuleExplicit
from pylith.utils.profiling import resourceUsageString

# Explicit class
class Explicit(Formulation, ModuleExplicit):
  """
  Python Explicit object for solving equations using an explicit
  formulation.

  The formulation has the general form, [A(t)] {u(t+dt)} = {b(t)},
  where we want to solve for {u(t+dt)}, A(t) is usually constant
  (i.e., independent of time), and {b(t)} usually depends on {u(t)}
  and {u(t-dt)}.

  Jacobian: A(t)
  solution: u(t+dt)
  residual: b(t) - A(t) \hat u(t+dt)
  constant: b(t)

  Factory: pde_formulation.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Formulation.Inventory):
    """
    Python object for managing Explicit facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Formulation facilities and properties.
    ##
    ## \b Properties
    ## @li \b norm_viscosity Normalized viscosity for numerical damping.
    ##
    ## \b Facilities
    ## @li \b solver Algebraic solver.

    import pyre.inventory

    normViscosity = pyre.inventory.float("norm_viscosity", default=0.1)
    normViscosity.meta['tip'] = "Normalized viscosity for numerical damping."

    from SolverLumped import SolverLumped
    solver = pyre.inventory.facility("solver", family="solver",
                                     factory=SolverLumped)
    solver.meta['tip'] = "Algebraic solver."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="explicit"):
    """
    Constructor.
    """
    Formulation.__init__(self, name)
    ModuleExplicit.__init__(self)
    self._loggingPrefix = "TSEx "
    self.dtStable = None
    return


  def elasticityIntegrator(self):
    """
    Get integrator for elastic material.
    """
    from pylith.feassemble.ElasticityExplicit import ElasticityExplicit
    integrator = ElasticityExplicit()
    integrator.normViscosity(self.normViscosity)
    return integrator


  def initialize(self, dimension, normalizer):
    """
    Initialize problem for explicit time integration.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)
    
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    self._initialize(dimension, normalizer)

    #from pylith.utils.petsc import MemoryLogger
    #memoryLogger = MemoryLogger.singleton()
    #memoryLogger.setDebug(0)
    #memoryLogger.stagePush("Problem")

    # Allocate other fields, reusing layout from solnIncr
    if 0 == comm.rank:
      self._info.log("Creating other fields.")
    self.fields.add("soln(t-dt)", "solution")
    self.fields.add("solnDeriv1(t)", "solution_deriv")
    self.fields.add("solnDeriv2(t)", "solution_deriv2")
    self.fields.copyLayout("solnIncr(t->t+dt)")
    self._debug.log(resourceUsageString())

    # Setup fields and set to zero
    dispTmdt = self.fields.get("soln(t-dt)")
    dispTmdt.zeroAll()
    dispT = self.fields.get("soln(t)")
    dispT.zeroAll()
    residual = self.fields.get("residual")
    residual.zeroAll()
    residual.createScatter(residual.mesh())

    lengthScale = normalizer.lengthScale()
    timeScale = normalizer.timeScale()
    velocityScale = lengthScale / timeScale
    velocityT = self.fields.get("velocity(t)")
    velocityT.scale(velocityScale.value)
    velocityT.zeroAll()

    accelerationScale = velocityScale / timeScale
    accelerationT = self.fields.get("acceleration(t)")
    accelerationT.scale(accelerationScale.value)
    accelerationT.zeroAll()

    self._debug.log(resourceUsageString())
    #memoryLogger.stagePop()

    if 0 == comm.rank:
      self._info.log("Creating lumped Jacobian matrix.")
    from pylith.topology.Field import Field
    jacobian = Field(self.mesh())
    jacobian.newSection(jacobian.VERTICES_FIELD, dimension)
    jacobian.allocate()
    jacobian.label("jacobian")
    jacobian.vectorFieldType(jacobian.VECTOR)
    self.jacobian = jacobian
    self._debug.log(resourceUsageString())

    #memoryLogger.stagePush("Problem")
    if 0 == comm.rank:
      self._info.log("Initializing solver.")
    self.solver.initialize(self.fields, self.jacobian, self)
    self._debug.log(resourceUsageString())

    #memoryLogger.stagePop()
    #memoryLogger.setDebug(0)
    self._eventLogger.eventEnd(logEvent)
    return


  def prestep(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    logEvent = "%sprestep" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)
    
    solnIncr = self.fields.get("solnIncr(t->t+dt)")
    for constraint in self.constraints:
      constraint.setFieldIncr(t, t+dt, solnIncr)

    needNewJacobian = False
    for integrator in self.integrators:
      integrator.timeStep(dt)
      if integrator.needNewJacobian():
        needNewJacobian = True
    if needNewJacobian:
      self._reformJacobian(t, dt)

    self._eventLogger.eventEnd(logEvent)
    return


  def step(self, t, dt):
    """
    Advance to next time step.
    """
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    self._reformResidual(t, dt)
    
    if 0 == comm.rank:
      self._info.log("Solving equations.")

    residual = self.fields.get("residual")
    solnIncr = self.fields.get("solnIncr(t->t+dt)")
    self.solver.solve(solnIncr, self.jacobian, residual)

    return


  def poststep(self, t, dt):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)
    
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    # The velocity and acceleration at time t depends on the
    # displacement at time t+dt, we want to output BEFORE updating the
    # displacement fields so that the displacement, velocity, and
    # acceleration files are all at time t.
    if 0 == comm.rank:
      self._info.log("Writing solution fields.")
    for output in self.output.components():
      output.writeData(t, self.fields)
    self._writeData(t)

    # Update displacement field from time t to time t+dt.
    solnIncr = self.fields.get("solnIncr(t->t+dt)")
    solnT = self.fields.get("soln(t)")
    solnTmdt = self.fields.get("soln(t-dt)")

    solnTmdt.copy(solnT)
    solnT.add(solnIncr)
    solnIncr.zeroAll()

    # Complete post-step processing.
    Formulation.poststep(self, t, dt)

    self._eventLogger.eventEnd(logEvent)    
    return


  def prestepElastic(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()
    
    if 0 == comm.rank:
      self._info.log("Setting constraints.")
    soln = self.fields.get("solnIncr(t->t+dt)")
    soln.zeroAll()
    for constraint in self.constraints:
      constraint.setField(t+dt, soln)

    needNewJacobian = False
    for integrator in self.integrators:
      integrator.timeStep(dt)
      if integrator.needNewJacobian():
        needNewJacobian = True
    if needNewJacobian:
      self._reformJacobian(t, dt)

    return


  def getTimeStep(self):
    """
    Get stable time step for advancing forward in time. Use cached
    value if available.

    Assume stable time step depends only on initial elastic properties
    and original mesh geometry.
    """
    logEvent = "%stimestep" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    if self.dtStable is None:
      self.dtStable = self.timeStep.timeStep(self.mesh(), self.integrators)
    self._eventLogger.eventEnd(logEvent)
    return self.dtStable
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Formulation._configure(self)

    self.normViscosity = self.inventory.normViscosity
    self.solver = self.inventory.solver
    return


  def _reformJacobian(self, t, dt):
    """
    Reform Jacobian matrix for operator.
    """
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    self._debug.log(resourceUsageString())
    if 0 == comm.rank:
      self._info.log("Integrating Jacobian operator.")
    self._eventLogger.stagePush("Reform Jacobian")

    self.updateSettings(self.jacobian, self.fields, t, dt)
    ModuleExplicit.reformJacobianLumped(self)

    self._eventLogger.stagePop()

    if self.viewJacobian:
      self.jacobian.view("Jacobian")

    self._debug.log(resourceUsageString())
    return


# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Explicit.
  """
  return Explicit()


# End of file 
