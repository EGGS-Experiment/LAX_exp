from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class Cleanup(LAXSubsequence):
    """
    Subsequence: Cleanup

    Reset device states to their default for client use after a sequence has been run.
    """
    name = 'cleanup'

    def build_subsequence(self):
        # get devices
        self.setattr_device('urukul0_cpld')
        self.setattr_device('urukul1_cpld')
        self.setattr_device('urukul2_cpld')
        self.setattr_device('ttl12')
        self.setattr_device('ttl13')
        self.setattr_device('ttl14')

    @kernel(flags={"fast-math"})
    def run(self):
        # reset core device, RTIOs, and FIFOs
        self.core.reset()

        # reset hardware to allow use
        with parallel:
            # enable all RF switches
            with parallel:
                self.ttl12.off()
                self.ttl13.off()
                self.ttl14.off()

            # reset qubit board
            with sequential:
                self.urukul0_cpld.set_profile(0)
                self.urukul0_cpld.io_update.pulse_mu(8)
                self.urukul0_cpld.cfg_switches(0b0000)

            # reset motional board to rescue parameters
            with sequential:
                self.urukul1_cpld.set_profile(0)
                self.urukul1_cpld.io_update.pulse_mu(8)
                self.urukul1_cpld.cfg_switches(0b0000)

            # reset main board to rescue parameters
            with sequential:
                self.urukul2_cpld.set_profile(0)
                self.urukul2_cpld.io_update.pulse_mu(8)
                self.urukul2_cpld.cfg_switches(0b1110)

        self.core.break_realtime()
