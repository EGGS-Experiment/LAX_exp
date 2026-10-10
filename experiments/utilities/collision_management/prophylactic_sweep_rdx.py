from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from numpy import int32
from LAX_exp.language import *
from LAX_exp.system.subsequences import RescueIon


class ProphylacticSweepRDX(LAXExperiment, Experiment):
    """
    Utility: Prophylactic Sweep RDX

    Apply a tickle on the RF rods as a post-exposure prophylaxis against
        undesired mass species from being trapped.
    Note: this utility does not set the 397nm beams, so user configs from EGGS GUI will
        persist during operation.
    """
    name = 'Prophylactic Sweep RDX'
    kernel_invariants = {
        "profile_tickle", "profile_parametric",
        "ampl_parametric_asf", "att_parametric_mu", "ampl_tickle_asf", "att_tickle_mu",

        "mod_freq_mu_list", "mod_time_mu", "time_cooling_holdoff_mu"
    }

    def build_experiment(self):
        # allocate relevant DDS profiles
        self.profile_tickle =      0
        self.profile_parametric =  0

        # modulation - general
        self.setattr_argument("repetitions",        NumberValue(default=5, precision=0, step=1, min=1, max=10000),
                              tooltip="Number of times to sweep through mod_freq_khz_list.")
        self.setattr_argument("mod_time_ms",        NumberValue(default=20, precision=3, step=1, min=0.001, max=10000000, scale=1., unit="ms"),
                              tooltip="Amount of time to apply prophylaxis at a given frequency.")
        self.setattr_argument("mod_freq_khz_list",  Scannable(
                                                            default=[
                                                                CenterScan(1686, 1000, 100., randomize=True),
                                                                RangeScan(1800, 2300, 501, randomize=True),
                                                            ],
                                                            global_min=1, global_max=100000, global_step=1.,
                                                            unit="kHz", scale=1, precision=3
                                                        ),
                              tooltip="List of frequencies to apply prophylaxis at."
                                      "These should simply be the target secular frequencies "
                                      "(since scaling can be separately applied to the tickle and parametric DDSs.")

        # modulation - tickle config
        self.setattr_argument("ampl_tickle_pct",    NumberValue(default=30, precision=3, step=5, min=0.01, max=100, scale=1., unit="%"),
                              group="{}.tickle".format("spegra"),
                              tooltip="DDS amplitude to use for tickling.")
        self.setattr_argument("att_tickle_db",      NumberValue(default=3, precision=1, step=0.5, min=0., max=31.5, scale=1., unit="dB"),
                              group="{}.tickle".format("spegra"),
                              tooltip="DDS attenuation to use for tickling.")
        self.setattr_argument("freq_tickle_scale",  NumberValue(default=1, precision=6, step=0.1, min=0., max=100., scale=1.),
                              group="{}.tickle".format("spegra"),
                              tooltip="Frequency scale for the tickle DDS.\n"
                                      "Scales all frequencies in mod_freq_khz_list by this value for the tickle DDS only.\n"
                                      "Useful when applying tickle AND parametric to get sharper excitation.")

        # modulation - parametric config
        self.setattr_argument("ampl_parametric_pct",    NumberValue(default=30, precision=3, step=5, min=0.01, max=100, scale=1., unit="%"),
                              group="{}.parametric".format("spegra"),
                              tooltip="DDS amplitude to use for parametric.")
        self.setattr_argument("att_parametric_db",      NumberValue(default=3, precision=1, step=0.5, min=0., max=31.5, scale=1., unit="dB"),
                              group="{}.parametric".format("spegra"),
                              tooltip="DDS attenuation to use for parametric.")
        self.setattr_argument("freq_parametric_scale",  NumberValue(default=2, precision=6, step=0.1, min=0., max=100., scale=1.),
                              group="{}.parametric".format("spegra"),
                              tooltip="Frequency scale for the parametric DDS.\n"
                                      "Scales all frequencies in mod_freq_khz_list by this value for the parametric DDS only.\n"
                                      "Useful when applying tickle AND parametric to get sharper excitation.")

        # get relevant devices
        self.setattr_device("dds_dipole")
        self.setattr_device("dds_parametric")

        # initialize other subsequences
        self.rescue_subsequence = RescueIon(self)

    def prepare_experiment(self):
        # convert DDS parameters to machine units
        self.ampl_tickle_asf =  self.dds_dipole.amplitude_to_asf(self.ampl_tickle_pct / 100.)
        self.att_tickle_mu =    self.dds_dipole.cpld.att_to_mu(self.att_tickle_db * dB)

        self.ampl_parametric_asf =  self.dds_parametric.amplitude_to_asf(self.ampl_parametric_pct / 100.)
        self.att_parametric_mu =    self.dds_parametric.cpld.att_to_mu(self.att_parametric_db * dB)

        # set up modulation timings etc.
        self.mod_freq_mu_list = [self.dds_parametric.frequency_to_ftw(freq_khz * kHz)
                                 for freq_khz in self.mod_freq_khz_list]
        self.mod_time_mu =      self.core.seconds_to_mu(self.mod_time_ms * ms)
        self.time_cooling_holdoff_mu = self.core.seconds_to_mu(3 * ms)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.mod_freq_mu_list),
                2)


    """
    MAIN SEQUENCE
    """
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        self.dds_dipole.set_att_mu(self.att_tickle_mu)
        self.dds_parametric.set_att_mu(self.att_parametric_mu)

        self.dds_dipole.set_profile(self.profile_tickle)
        self.dds_parametric.set_profile(self.profile_parametric)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        for i in range(self.repetitions):
            for j in range(len(self.mod_freq_mu_list)):
                # get target freq
                freq_mu = self.mod_freq_mu_list[j]

                # add holdoff for recooling the ion and synchronize with timeline
                delay_mu(self.time_cooling_holdoff_mu)
                self.core.wait_until_mu(now_mu())

                # prepare DDS waveforms
                self.core.break_realtime()
                self.dds_dipole.set_mu(int32(freq_mu * self.freq_tickle_scale),
                                       asf=self.ampl_tickle_asf,
                                       profile=self.profile_tickle,
                                       phase_mode=PHASE_MODE_CONTINUOUS)
                self.dds_parametric.set_mu(int32(freq_mu * self.freq_parametric_scale),
                                           asf=self.ampl_parametric_asf,
                                           profile=self.profile_parametric,
                                           phase_mode=PHASE_MODE_CONTINUOUS)

                # run excitation (tickle + parametric) for target time
                self.dds_dipole.on()
                self.dds_parametric.on()
                delay_mu(self.mod_time_mu)
                self.dds_dipole.off()
                self.dds_parametric.off()

                # clean up & periodically check termination
                self.update_results(freq_mu, freq_mu)
                if j % 50 == 0: self.check_termination()

