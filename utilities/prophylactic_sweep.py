import numpy as np
from artiq.experiment import *


class ProphylacticSweep(EnvExperiment):
    """
    Prophylactic Sweep

    Apply a tickle on top of the RF rods as a prophylaxis against undesired
    mass species from being trapped.
    """


    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # modulation
        self.setattr_argument("mod_time_s",                         NumberValue(default=10, ndecimals=3, step=1, min=0.001, max=10000000))
        self.setattr_argument("mod_att_db",                         NumberValue(default=31, ndecimals=1, step=0.5, min=0, max=31.5))
        self.setattr_argument("mod_freq_mhz_list",                  Scannable(
                                                                        default=CenterScan(1.686, 0.01, 0.0001, randomize=True),
                                                                        global_min=0, global_max=1000, global_step=0.001,
                                                                        unit="MHz", scale=1, ndecimals=5
                                                                    ))

    def prepare(self):
        # modulation control
        self.mod_dds =                                              self.get_device("urukul1_ch2")
        self.mod_dds_ampl_pct =                                     self.mod_dds.amplitude_to_asf(0.35)
        self.mod_dds_att_mu =                                       self.mod_dds.cpld.att_to_mu(self.mod_att_db * dB)
        self.mod_freq_mu_list =                                     np.array([
                                                                        self.mod_dds.frequency_to_ftw(freq_mhz * MHz)
                                                                        for freq_mhz in self.mod_freq_mhz_list
                                                                    ])
        self.mod_time_mu =                                          self.core.seconds_to_mu(self.mod_time_s)

        # cooling holdoff time
        self.time_cooling_holdoff_mu =                              self.core.seconds_to_mu(3 * ms)


    @kernel(flags={"fast-math"})
    def run(self):
        # reset core device
        self.core.reset()

        # prepare devices for experiment
        self.prepareDevices()
        self.core.break_realtime()


        # MAIN LOOP
        # sweep modulation frequency
        for freq_mu in self.mod_freq_mu_list:

            # add holdoff period for recooling the ion
            delay_mu(self.time_cooling_holdoff_mu)

            # set modulation frequency and tickle
            self.mod_dds.set_mu(freq_mu, asf=self.mod_dds_ampl_pct)
            self.mod_dds.cfg_sw(True)

            # apply prophylaxis for given time
            delay_mu(self.mod_time_mu)

            # stop counting and upload
            self.mod_dds.cfg_sw(False)

            # add delay at end
            self.core.break_realtime()


    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        # configure rf modulation source
        self.mod_dds.cfg_sw(False)
        self.mod_dds.set_att_mu(self.mod_dds_att_mu)
