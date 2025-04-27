import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class CharacteristicReadout(LAXSubsequence):
    """
    Subsequence: Characteristic Readout

    Directly reconstruct the Characteristic function of a given motional state using the Fluhmann technique
    (https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.043602).
    """
    name = 'characteristic_readout'
    kernel_invariants = {
        # core hardware
        'singlepass0', 'singlepass1',

        # hardware parameters
        'freq_singlepass_default_ftw_list', 'ampl_singlepass_default_asf_list', 'att_singlepass_default_mu_list',
        'ampl_doublepass_default_asf', 'att_doublepass_default_mu',
        'freq_sigmax_ftw', 'ampl_sigmax_asf', 'att_sigmax_mu', 'time_sigmax_mu',
        'time_force_herald_slack_mu', 'profile_729_target',

        # cat-state parameters
        'ampls_cat_asf', 'atts_cat_mu', 'time_pulse1_cat_mu', 'phases_pulse1_cat_pow', 'phase_characteristic_axis_pow',
        'phases_pulse5_cat_pow', 'phases_pulse5_cat_update_dir'
    }

    def build_subsequence(self, profile_dds: TInt32 = 0) -> TNone:
        """
        Defines the main interface for the subsequence.
        Arguments:
            fdfd
        """
        # get relevant devices
        self.setattr_device('qubit')
        self.setattr_device('repump_qubit')
        self.setattr_device('pump')

        # set subsequence parameters
        self.profile_dds = profile_dds

        # configure detection of real/imag part of characteristic function
        self.setattr_argument("characteristic_axis", EnumerationValue(['Both', 'Real', 'Imaginary'], default='Both'),
                              group=self.name,
                              tooltip="Selects the real/imag component of the characteristic function by either applying a sigma_x operation (Imag), or not (Real)."
                                      "The 'Both' option enables measurement of both real and imag components within a single experiment.")
        self.setattr_argument("phase_characteristic_axis_turns", NumberValue(default=0., precision=3, step=0.1, min=-1.0, max=1.0),
                              group=self.name,
                              tooltip="Sets the relative phase of the sigma_x operation used to define the real/imag axis of the characteristic function.")

        # pulse 5 - cat 1: characteristic readout protocol
        self.setattr_argument("target_char_read_phase", EnumerationValue(['RSB', 'BSB', 'RSB-BSB', 'RSB+BSB'], default='RSB-BSB'),
                              group=self.name,
                              tooltip="Selects how phase is controlled when scanning over the plane of the"
                                      "characteristic function.")
        self.setattr_argument("phases_char_read_turns", PYONValue([0., 0.]), group=self.name,
                              tooltip="[rsb_turns (ch1), bsb_turns (ch2)]")
        self.setattr_argument("time_char_read_x_us_list", Scannable(
                                                            default=[
                                                                RangeScan(-500, 500, 10, randomize=True),
                                                                ExplicitScan([100]),
                                                            ],
                                                            global_min=-100000, global_max=100000, global_step=1,
                                                            unit="us", scale=1, precision=5), group=self.name)
        self.setattr_argument("time_char_read_y_us_list", Scannable(
                                                            default=[
                                                                RangeScan(-500, 500, 10, randomize=True),
                                                                ExplicitScan([100]),
                                                            ],
                                                            global_min=-100000, global_max=100000, global_step=1,
                                                            unit="us", scale=1, precision=5), group=self.name)

    def prepare_subsequence(self):
        """
        Prepare & precompute experimental values.
        """
        '''
        CONVERT VALUES TO MACHINE UNITS - DEFAULTS
        '''
        # defaults - singlepass AOM
        self.singlepass0 = self.get_device("urukul0_ch1")
        self.singlepass1 = self.get_device("urukul0_ch2")
        self.freq_singlepass_default_ftw_list = [self.singlepass0.frequency_to_ftw(freq_mhz * MHz)
                                                 for freq_mhz in self.freq_singlepass_default_mhz_list]
        self.ampl_singlepass_default_asf_list = [self.singlepass0.amplitude_to_asf(ampl_asf / 100.)
                                                 for ampl_asf in self.ampl_singlepass_default_pct_list]
        self.att_singlepass_default_mu_list =   [att_to_mu(att_db * dB)
                                                 for att_db in self.att_singlepass_default_db_list]

        # defaults - doublepass AOM
        self.ampl_doublepass_default_asf =     self.qubit.amplitude_to_asf(self.ampl_doublepass_default_pct / 100.)
        self.att_doublepass_default_mu =       att_to_mu(self.att_doublepass_default_db * dB)

        # defaults - sigma_x pulses
        self.freq_sigmax_ftw =  self.qubit.frequency_to_ftw(self.freq_sigmax_mhz * MHz)
        self.ampl_sigmax_asf =  self.qubit.amplitude_to_asf(self.ampl_sigmax_pct / 100.)
        self.att_sigmax_mu =    att_to_mu(self.att_sigmax_db * dB)
        self.time_sigmax_mu =   self.core.seconds_to_mu(self.time_sigmax_us * us)

        # defaults - cat
        freq_cat_center_ftw_list =  np.array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                              for freq_mhz in self.freq_cat_center_mhz_list])
        freq_cat_secular_ftw_list = np.array([self.singlepass0.frequency_to_ftw(freq_khz * kHz)
                                              for freq_khz in self.freq_cat_secular_khz_list])
        self.ampls_cat_asf =    np.array([self.singlepass0.amplitude_to_asf(ampl_pct / 100.)
                                          for ampl_pct in self.ampls_cat_pct], dtype=np.int32)
        self.atts_cat_mu =      np.array([att_to_mu(att_db * dB)
                                          for att_db in self.atts_cat_db], dtype=np.int32)

        '''
        CONVERT VALUES TO MACHINE UNITS - PULSES
        '''
        # define characteristic axis (via sigma_x)
        self.phase_characteristic_axis_pow =  self.qubit.turns_to_pow(self.phase_characteristic_axis_turns)
        # configure whether real/imag/both axes of the characteristic function are to be measured
        if self.characteristic_axis == "Real":          characteristic_axis_list = [False]
        elif self.characteristic_axis == "Imaginary":   characteristic_axis_list = [True]
        elif self.characteristic_axis == "Both":        characteristic_axis_list = [True, False]

        # pulse 5 - cat 1
        self.phases_pulse5_cat_pow = np.array([self.singlepass0.turns_to_pow(phas_pow)
                                               for phas_pow in self.phases_pulse5_cat_turns], dtype=np.int32)

        if self.target_pulse5_cat_phase == 'RSB':
            self.phases_pulse5_cat_update_dir = np.array([1, 0], dtype=np.int32)
        elif self.target_pulse5_cat_phase == 'BSB':
            self.phases_pulse5_cat_update_dir = np.array([0, 1], dtype=np.int32)
        elif self.target_pulse5_cat_phase == 'RSB-BSB':
            self.phases_pulse5_cat_update_dir = np.array([1, -1], dtype=np.int32)
        elif self.target_pulse5_cat_phase == 'RSB+BSB':
            self.phases_pulse5_cat_update_dir = np.array([1, 1], dtype=np.int32)

        # create sampling grid in radial coordinates
        vals_pulse5_mu_pow_list = np.array([
            [
                self.core.seconds_to_mu(math.sqrt(x_us ** 2. + y_us ** 2.) * us),
                self.singlepass0.turns_to_pow(np.arctan2(y_us, x_us) / (2. * np.pi))
            ]
            for x_us in self.time_pulse5_cat_x_us_list
            for y_us in self.time_pulse5_cat_y_us_list
        ], dtype=np.int64)

        '''
        CREATE EXPERIMENT CONFIG
        '''
        # use generator to flatten list with some tuples
        def flatten(xs):
            for x in xs:
                if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
                    yield from flatten(x)
                else:
                    yield x

        # create an array of values for the experiment to sweep
        self.config_experiment_list = np.array([
            list(flatten(vals))
            for vals in product(
                freq_cat_center_ftw_list, freq_cat_secular_ftw_list,
                vals_pulse5_mu_pow_list,
                characteristic_axis_list
            )
        ], dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        pass

    @kernel(flags={"fast-math"})
    def run_time(self, time_char_read_mu: TInt64, phase_char_read_pow: TInt32, char_ax: TBool,
                 time_start_mu: TInt64) -> TNone:
        """
        TODO: DOCUMENT
        :param time_readout_mu:
        :return:
        """
        # set variables for
        char_read_phases_pow = [
            self.phases_pulse5_cat_pow[0] + self.phases_pulse5_cat_update_dir[0] * phase_char_read_pow,
            self.phases_pulse5_cat_pow[1] + self.phases_pulse5_cat_update_dir[1] * phase_char_read_pow,
        ]

        # prepare spin state for characteristic readout
        # note: need to set correct profile for normal quenching
        # otherwise might be stuck in SBC quench params)
        self.pump.readout()
        self.repump_qubit.on()
        delay_mu(self.initialize_subsequence.time_repump_qubit_mu)
        self.repump_qubit.off()

        # select real vs imag part of characteristic function (via a pi/2 pulse)
        if char_ax:
            self.pulse_sigmax(time_start_mu, self.phase_characteristic_axis_pow)

        # pulse 5: cat #2
        self.pulse_bichromatic(time_start_mu, time_char_read_mu,
                               char_read_phases_pow,
                               freq_cat_center_ftw, freq_cat_secular_ftw)

    '''
    HELPER FUNCTIONS
    '''
    @kernel(flags={"fast-math"})
    def pulse_bichromatic(self, time_start_mu: TInt64, time_pulse_mu: TInt64, phas_pow_list: TList(TInt32),
                          freq_carrier_ftw: TInt32, freq_secular_ftw: TInt32) -> TNone:
        """
        Run a phase-coherent bichromatic pulse on the qubit.
        Arguments:
            time_start_mu: fiducial timestamp for initial start reference (in machine units).
            time_pulse_mu: length of pulse (in machine units).
            phas_pow_list: relative phase offset for the beams (RSB, BSB) (in pow).
            freq_carrier_ftw: carrier frequency (set by the double pass) in FTW.
            freq_secular_ftw: bichromatic separation frequency (from central frequency) in FTW.
        """
        # set up relevant beam waveforms
        self.qubit.set_mu(
            freq_carrier_ftw, asf=self.ampl_sigmax_asf, pow_=0,
            profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.singlepass0.set_mu(
            self.freq_singlepass_default_ftw_list[0]-freq_secular_ftw, asf=self.ampls_cat_asf[0],
            pow_=phas_pow_list[0], profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.singlepass1.set_mu(
            self.freq_singlepass_default_ftw_list[1]+freq_secular_ftw, asf=self.ampls_cat_asf[1],
            pow_=phas_pow_list[1], profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )

        # set all attenuators together
        a = self.qubit.cpld.att_reg & ~(
                (0xFF << (0 * 8)) |
                (0xFF << (1 * 8)) |
                (0xFF << (2 * 8))
        )
        a |= (
                (self.att_doublepass_default_mu << (0 * 8)) |
                (self.atts_cat_mu[0] << (1 * 8)) |
                (self.atts_cat_mu[1] << (2 * 8))
        )
        self.qubit.cpld.set_all_att_mu(a)

        # run bichromatic pulse
        self.singlepass0.sw.on()
        self.singlepass1.sw.on()
        self.qubit.on()
        delay_mu(time_pulse_mu)
        self.qubit.off()
        self.singlepass1.sw.off()

