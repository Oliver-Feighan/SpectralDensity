from openmm.app import DCDFile
from openmm.unit import nanometer

import mdtraj

class CHLReporter(object):
    """DCDReporter outputs a series of frames from a Simulation to a DCD file.
    To use it, create a DCDReporter, then add it to the Simulation's list of reporters.
    """


    def __init__(self, file, reportInterval, append=False, enforcePeriodicBox=None):
        """Create a DCDReporter.
        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        append : bool=False
            If True, open an existing DCD file to append to.  If False, create a new file.
        enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
        """
        self._reportInterval = reportInterval
        self._append = append
        self._enforcePeriodicBox = enforcePeriodicBox
        if append:
            mode = 'r+b'
        else:
            mode = 'wb'
        self._out = open(file, mode)
        self._dcd = None


    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.
        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        Returns
        -------
        tuple
            A six element tuple. The first element is the number of steps
            until the next report. The next four elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.  The final element specifies whether
            positions should be wrapped to lie in a single periodic box.
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, True, False, False, False, self._enforcePeriodicBox)


    def report(self, simulation, state):
        """Generate a report.
        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """

        if self._dcd is None:
            
            full_top = mdtraj.Topology.from_openmm(simulation.topology)
            self._chl_indices = full_top.select("resname =~ 'BCL*'")
            chl_top = full_top.subset(self._chl_indices) 
            openmm_chl_top = chl_top.to_openmm()
            self._dcd = DCDFile(
                self._out, openmm_chl_top, simulation.integrator.getStepSize(),
                simulation.currentStep, self._reportInterval, self._append
            )
        crds = state.getPositions(asNumpy=True)
        crds = crds[self._chl_indices]
        self._dcd.writeModel(crds, periodicBoxVectors=state.getPeriodicBoxVectors())

    def __del__(self):
        self._out.close()
