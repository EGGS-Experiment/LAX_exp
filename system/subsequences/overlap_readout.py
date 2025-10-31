from numpy import int64
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS, RAM_MODE_RAMPUP

from LAX_exp.language import *
from LAX_exp.system.subsequences import QubitRAP


class OverlapReadout(QubitRAP):
    """
    Subsequence: Overlap Readout

    Inherits from QubitRAP class.
    Use the motional overlap technique from F.Wolf (P.O. Schmidt group) to directly interrogate
        the population of a given fock state (https://www.nature.com/articles/s41467-019-10576-4).
    The maximum interrogable fock state is determined by the number of available auxiliary states.
    """
    name = 'overlap_readout'
    kernel_invariants = {
        # configuration
        "profile_shelve", "num_pulses",

        # hardware values
        "att_overlap_mu", "freq_carrier_ftw", "ampl_carrier_asf", "time_carrier_mu",
        "freq_rap_center_ftw", "freq_rap_dev_ftw", "time_rap_mu"
    }

    def build_subsequence(self, ram_profile: TInt32 = -1, profile_shelve: TInt32 = -1,
                          ram_addr_start: TInt32 = 0x00, num_samples: TInt32 = 200,
                          pulse_shape: TStr = "blackman"):
        """
        Defines the main interface for the subsequence.
        :param ram_profile: the AD9910 RAM profile to use for RAP + pulse shaping.
        :param profile_shelve: the AD9910 profile (non-RAM) to use for shelving.
        :param ram_addr_start: the beginning RAM register address for pulse shaping.
            Must be in [0, 923].
        :param num_samples: the number of samples to use for the pulse shape.
            Must result in a final RAM address <= 1023.
        :param pulse_shape: the pulse shape to use. Must be supported by available_pulse_shapes.
        """
        # extend our kernel_invariants with parent's kernel_invariants
        # (since we redefine them here)
        kernel_invariants_parent = getattr(super(), "kernel_invariants", set())
        self.kernel_invariants = self.kernel_invariants | kernel_invariants_parent

        # extract desired args before passing them onto parent (i.e. QubitRAP)
        self.profile_shelve = profile_shelve
        _argstr = "overlap" # create short string for argument grouping

        # overlap readout: arguments
        self.setattr_argument("fock_state_readout", EnumerationValue(["0", "1", "2", "3"], default="0"), group="{}".format(_argstr))
        self.setattr_argument("shelving_config",    PYONValue({0: (113.3558, 50., 5.08), 1: (104.1466, 50., 12.26), 2: (110.2858, 50., 25.40)}),
                              group="{}".format(_argstr), tooltip="{fock_state_n: (freq_mhz, ampl_pct, time_us).")
        self.setattr_argument("att_overlap_db",     NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, unit="dB", scale=1.),
                              group="{}".format(_argstr),
                              tooltip="Same attenuation will be used for all overlap readout -"
                                      "i.e. RAP, shelving, and carrier pulses.")

        # carrier: arguments
        self.setattr_argument("freq_carrier_mhz",   NumberValue(default=101.0781, precision=6, step=1e-2, min=60, max=200, unit="MHz", scale=1.),
                              group="{}.carr".format(_argstr))
        self.setattr_argument("ampl_carrier_pct",   NumberValue(default=50., precision=3, step=5, min=1, max=50, unit="%", scale=1.),
                              group="{}.carr".format(_argstr))
        self.setattr_argument("time_carrier_us",    NumberValue(default=2.23, precision=3, min=1, max=1e5, step=1, unit="us", scale=1.),
                              group="{}.carr".format(_argstr))

        # RAP: arguments
        self.setattr_argument("freq_rap_center_mhz",    NumberValue(default=100.7471, precision=6, step=1e-2, min=60, max=200, unit="MHz", scale=1.),
                              group="{}.RAP".format(_argstr))
        self.setattr_argument("freq_rap_dev_khz",       NumberValue(default=100., precision=2, step=0.01, min=1, max=1e4, unit="kHz", scale=1.),
                              group="{}.RAP".format(_argstr),
                              tooltip="Single-sided frequency deviation for RAP."
                                      "Chirp will be done from [freq_center - freq_rap_dev, freq_center + freq_rap_dev].")
        self.setattr_argument("time_rap_us",    NumberValue(default=500., precision=3, min=1, max=1e5, step=1, unit="us", scale=1.),
                              group="{}.RAP".format(_argstr))
        self.setattr_argument("ampl_rap_pct",   NumberValue(default=50., precision=3, step=5, min=1, max=50, unit="%", scale=1.),
                              group="{}.RAP".format(_argstr))

        # build parent subsequence (QubitRAP)
        super().build_subsequence(
            ram_profile=ram_profile, ram_addr_start=ram_addr_start, num_samples=num_samples,
            ampl_max_pct=self.ampl_rap_pct, pulse_shape=pulse_shape
        )

    def prepare_subsequence(self):
        """
        Prepare and precompute experiment values for speedy evaluation.
        """
        # call parent prepare_subsequence
        super().prepare_subsequence()
        # question: if I do super().prepare_subsequence, does it call their _prepare_argument_checks,
        #  or MY _prepare_argument_checks?
        # answer: it runs MY _prepare_argument_checks (reasonably)

        '''CONVERT PULSE CONFIGURATION'''
        # process fock state
        # note: this is OK b/c EnumerationValue limits possible inputs, and also b/c we check validity
        # in _prepare_argument_checks
        self.num_pulses = int(self.fock_state_readout)

        # note: create separate arrays for freq/ampl & time to avoid int32 conversion later on
        # (b/c all variables have to have same type)
        # create dummy configs if not needed (to avoid artiq typing issues)
        if self.num_pulses == 0:
            self.dds_configs = [
                [self.qubit.frequency_to_ftw(config[0] * MHz),
                 self.qubit.amplitude_to_asf(config[1] / 100.)]
                for config in [[200., 10., 10.]]
            ]
            self.pulse_times_mu = [
                self.core.seconds_to_mu(config[2] * us)
                for config in [[200., 10., 10.]]
            ]
        else:
            self.dds_configs = [
                [self.qubit.frequency_to_ftw(config[0] * MHz),
                 self.qubit.amplitude_to_asf(config[1] / 100.)]
                for config in self.shelving_config.values()
            ]
            self.pulse_times_mu = [
                self.core.seconds_to_mu(config[2] * us)
                for config in self.shelving_config.values()
            ]

        '''CONVERT VALUES TO MACHINE UNITS'''
        # overlap
        self.att_overlap_mu = att_to_mu(self.att_overlap_db)

        # carrier
        self.freq_carrier_ftw = self.qubit.frequency_to_ftw(self.freq_carrier_mhz * MHz)
        self.ampl_carrier_asf = self.qubit.amplitude_to_asf(self.ampl_carrier_pct / 100.)
        self.time_carrier_mu =  self.core.seconds_to_mu(self.time_carrier_us * us)

        # RAP
        self.freq_rap_center_ftw = self.qubit.frequency_to_ftw(self.freq_rap_center_mhz * MHz)
        self.freq_rap_dev_ftw = self.qubit.frequency_to_ftw(self.freq_rap_dev_khz * kHz)
        self.time_rap_mu = self.core.seconds_to_mu(self.time_rap_us * us)
        self._time_rap_actual_mu = int64(0) # holder for RAP actual time (which we receive in initialize)

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        super()._prepare_argument_checks()

        # ensure fock_state_readout is valid
        if self.fock_state_readout not in ("0", "1", "2", "3"):
            raise ValueError("Invalid fock_state_readout: {:}. Must be in ['0', '1', '2', '3'].".format(self.fock_state_readout))

        # ensure target AD9910 profile is valid
        if self.profile_shelve not in range(0, 7):
            raise ValueError("Invalid AD9910 profile for overlap_readout: {:d}. Must be in [0, 7].".format(self.profile_shelve))

        # only require further shelving_config checks if we actually need them for higher |n> detection
        if self.fock_state_readout != "0":
            # ensure shelving config is valid
            if not isinstance(self.shelving_config, dict):
                raise ValueError("Invalid shelving config type: {:}. Must be of type dict.".format(self.shelving_config))

            # parse and check shelving config
            shelving_keys = self.shelving_config.keys()
            shelving_values = self.shelving_config.values()

            # checking config keys
            if not all((isinstance(k, int) and 0 <= k <= 2 for k in shelving_keys)):
                raise ValueError("Invalid shelving config keys: {:}. Must be dict with keys as ints.".format(shelving_keys))
            # ensure shelving config keys are correctly numbered
            elif not ((max(shelving_keys) + 1 == len(shelving_keys)) and (len(set(shelving_keys)) == len(shelving_keys))):
                raise ValueError("Invalid shelving config keys: {:}. Must be contiguous ints in [0, 2].".format(shelving_keys))

            # check config values are of correct type
            if not all((isinstance(config, tuple) and len(config) == 3 and
                        all((isinstance(val, (int, float)) for val in config))
                        for config in shelving_values)):
                raise ValueError("Invalid shelving config: {:}. Must be length-3 tuples of floats.".format(shelving_keys))
            # check configs all fall within valid range
            elif not all((2. <= config[0] <= 400. and 0.001 <= config[1] <= 50. and 0.04 <= config[2] <= 1000000.
                        for config in shelving_values)):
                raise ValueError("Invalid shelving config values: {:}. Values out of range.".format(shelving_values))

            # ensure enough shelving for fock state readout
            if max(shelving_keys) + 1 < int(self.fock_state_readout):
                raise ValueError("Insufficient shelving configs ({:}) for target fock_state_readout ({:}).".format(
                    max(shelving_keys), self.fock_state_readout))


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_subsequence(self) -> TNone:
        """
        Prepare the subsequence immediately before run.
        note: we copy over parent QubitRAP's initialize_subsequence b/c
            we need to also define our own, and ARTIQ doesn't support super()
            on coredevice.
        """
        '''FROM PARENT (QubitRAP)'''
        # disable RAM + DRG and set matched latencies
        self.qubit.set_cfr1(ram_enable=0)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.qubit.set_cfr2(matched_latency_enable=1)
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

        # write waveform to RAM via RAMWriter
        self.ram_writer.write(self.ampl_asf_pulseshape_list, self.ram_addr_start)
        delay_mu(25000)

        # set up RAM profile correctly after waveform uploaded
        self.qubit.set_profile_ram(
            start=self.ram_addr_start, end=self.ram_addr_stop,
            step=0xFFF,  # note: step size irrelevant since it's set in configure()
            profile=self.ram_profile, mode=RAM_MODE_RAMPUP
        )
        self.qubit.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        '''OUR ADDITIONS'''
        # configure RAP pulse ahead of time (b/c we never change it)
        # note: we store _time_rap_actual_mu b/c it gives us the REAL timing (the RAP pulse has large step sizes)
        self._time_rap_actual_mu = self.configure(self.time_rap_mu, self.freq_rap_center_ftw, self.freq_rap_dev_ftw)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        """
        Run complete motional overlap readout sequence.
        For best performance, should be recorded onto DMA.
        """
        # prepare beam (do here since RAP and shelving sequence all use same attenuation)
        self.qubit.set_att_mu(self.att_overlap_mu)

        # run first RAP (always necessary)
        self.run_rap(self._time_rap_actual_mu)

        # run multi-state shelving
        for i in range(self.num_pulses):
            self._shelving_pulse(i)

    @kernel(flags={"fast-math"})
    def _shelving_pulse(self, shelve_num: TInt32 = 0) -> TNone:
        """
        Run a shelving pulse to isolate the nth fock state population.
        :param shelve_num: the auxiliary state to use for shelving. Must be one of [0, 1, 2].
        """
        # prepare beam profile (need to do this here since we do RAP on a different profile)
        self.qubit.set_profile(self.profile_shelve)

        # run shelving pulse (|S,n> => |aux>)
        self.qubit.set_mu(self.dds_configs[shelve_num][0], asf=self.dds_configs[shelve_num][1],
                          profile=self.profile_shelve, phase_mode=PHASE_MODE_CONTINUOUS)
        self.qubit.on()
        delay_mu(self.pulse_times_mu[shelve_num])
        self.qubit.off()

        # run carrier pulse (|D,n+1> => |S,n+1>)
        self.qubit.set_mu(self.freq_carrier_ftw, asf=self.ampl_carrier_asf,
                          profile=self.profile_shelve, phase_mode=PHASE_MODE_CONTINUOUS)
        self.qubit.on()
        delay_mu(self.time_carrier_mu)
        self.qubit.off()

        # run RAP to isolate fock state (|S,n+1> => |D,n>)
        # note: (RAP takes care of profiles by itself)
        self.run_rap(self._time_rap_actual_mu)

