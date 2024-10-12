import numpy as np
from artiq.experiment import *


class ProphylacticSweep(EnvExperiment):
    """
    Utility: Prophylactic Sweep

    Apply a tickle on top of the RF rods as a prophylaxis against undesired
    mass species from being trapped.
    """


    def build(self):
        # get core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("scheduler")

        # modulation
        self.setattr_argument("repetitions",            NumberValue(default=5, precision=0, step=1, min=1, max=10000))
        self.setattr_argument("mod_time_total_s",       NumberValue(default=20, precision=3, step=1, min=0.001, max=10000000))
        self.setattr_argument("mod_att_db",             NumberValue(default=3, precision=1, step=0.5, min=0., max=31.5))
        self.setattr_argument("mod_freq_khz_list",      Scannable(
                                                            default=[
                                                                CenterScan(1686, 1000, 100., randomize=True),
                                                                RangeScan(1800, 2300, 501, randomize=True),
                                                            ],
                                                            global_min=1, global_max=100000, global_step=1.,
                                                            unit="kHz", scale=1, precision=3
                                                        ))

    def prepare(self):
        # modulation devices
        self.mod_dds =                      self.get_device("urukul1_ch1")
        self.mod_dds_sw =                   self.get_device("ttl11")
        self.mod_dds_ampl_pct =             self.mod_dds.amplitude_to_asf(0.35)
        self.mod_dds_att_mu =               self.mod_dds.cpld.att_to_mu(self.mod_att_db * dB)

        # timing/iteration
        self.mod_freq_mu_list =             np.array([
                                                self.mod_dds.frequency_to_ftw(freq_khz * kHz)
                                                for freq_khz in self.mod_freq_khz_list
                                            ])
        self.mod_time_mu =                  self.core.seconds_to_mu(self.mod_time_total_s / (len(self.mod_freq_mu_list) * self.repetitions))
        self.time_cooling_holdoff_mu =      self.core.seconds_to_mu(3 * ms)

        # for completion manager
        self._num_loops =                   len(self.mod_freq_mu_list) * self.repetitions


    @kernel(flags={"fast-math"})
    def run(self):
        # reset core device
        self.core.wait_until_mu(now_mu())
        self.core.reset()

        # prepare devices for experiment
        self.prepareDevices()
        _loop_counter = 0


        # MAIN LOOP
        for i in range(self.repetitions):
            # sweep modulation frequency
            for freq_mu in self.mod_freq_mu_list:

                # synchronize with timeline
                self.core.wait_until_mu(now_mu())
                self.core.break_realtime()

                # add holdoff period for recooling the ion
                delay_mu(self.time_cooling_holdoff_mu)

                # set modulation frequency and tickle
                self.mod_dds.set_mu(freq_mu, asf=self.mod_dds_ampl_pct)
                self.mod_dds.sw.on()
                self.mod_dds_sw.on()

                # apply prophylaxis for given time
                delay_mu(self.mod_time_mu)

                # turn off output
                self.mod_dds.sw.off()

                # update progress bar
                self.updateCompletion(_loop_counter)
                _loop_counter += 1

                # check termination
                if self.scheduler.check_termination():
                    self.cleanupDevices()
                    return

        # clean up devices
        self.cleanupDevices()


    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        # configure rf modulation source
        self.mod_dds.sw.off()
        self.mod_dds.set_att_mu(self.mod_dds_att_mu)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def cleanupDevices(self):
        """
        Clean up devices.
        """
        # delay experiment cancellation until sequence has finished
        self.core.break_realtime()

        # clean up DDS to prevent leakage
        self.mod_dds.set_att(31.5 * dB)
        self.mod_dds.set_mu(0x01, asf=0x01)
        # switch off all switches
        self.mod_dds_sw.off()
        self.mod_dds.cfg_sw(False)
        self.mod_dds.sw.off()
        self.core.break_realtime()

    @rpc(flags={"async"})
    def updateCompletion(self, i: TInt32) -> TNone:
        self.set_dataset('management.dynamic.completion_pct', round(i / self._num_loops * 100., 3),
                         broadcast=True, persist=True, archive=False)
