from artiq.experiment import *
from artiq.coredevice import ad9910

from numpy import array, int32, int64, zeros

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, Readout, ReadoutAdaptive,
    RescueIon, QubitRAP
)

from LAX_exp.system.objects.PulseShaper import available_pulse_shapes
from LAX_exp.system.objects.dds_pulse_shaper import DDSPulseShaper


class SecFreqDriftTracker(LAXExperiment, Experiment):
    """
    Experiment: SecFreqDriftTracker

    Create and characterize cat states with projective state preparation.
    """
    name = 'Secular Frequency Drift Tracker'
    kernel_invariants = {
        # subsequences & objects
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence',
        'rescue_subsequence', 'dds_pulse_shaper',

        # hardware values - tickle

        # hardware values - intensity servo
        'enable_servo_relock', 'time_servo_relock_mu',

        # configs
        'profile_729_SBC',
        'profile_tickle_RAM',
        'config_experiment_list',
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=50, precision=0, step=1, min=1, max=100000))

        # allocate relevant beam profiles
        self.profile_729_SBC = 0

        # allocate profiles for dds tickle
        self.profile_tickle_RAM = 0

        # get subsequences
        self.sidebandcool_subsequence = SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )
        self.initialize_subsequence = InitializeQubit(self)
        self.readout_subsequence = Readout(self)
        self.rescue_subsequence = RescueIon(self)

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')
        self.setattr_device('dds_dipole')
        self.setattr_device('ttl8')

        # set build arguments
        self._build_arguments_ion_parameters()
        self._build_arguments_tickle()
        self._build_arguments_rap()
        self._build_arguments_intensity_servo()

        self.rap_subsequences = []
        self.att_reg_readout_rap_list = []

    def _build_arguments_ion_parameters(self):
        """
        Build arguments for ion frequencies
        """
        _argstr = "ion_parameters"

        self.setattr_argument("freq_secular_khz_list",
                             PYONValue([710]),
                             group=_argstr,
                             tooltip="Secular frequency used for tickle pulse (in kHz) applied via the urukul dds.")

        self.setattr_argument("freq_carrier_mhz", NumberValue(
            default=100.481,
            min=60., max=400, step=1,
            unit="MHz", scale=1, precision=6
        ), group=_argstr,
                              tooltip="Carrier frequency of the ion.\n"
                                      "Note: this is applied via the main doublepass DDS.\n")

    def _build_arguments_rap(self):
        _argstr = 'rap'

        self.setattr_argument("att_rap_db_list",
                             PYONValue([8]),
                             group=_argstr,
                             tooltip="rap att.")

        self.setattr_argument("freq_rap_dev_khz_list",
                             PYONValue([72]),
                             group=_argstr,
                             tooltip="rap att.")

        self.setattr_argument("time_rap_us_list",
                             PYONValue([1000]),
                             group=_argstr,
                             tooltip="rap att.")

    def _build_arguments_tickle(self):
        """
        Build core sweep arguments for the tickle pulse.
        """
        _argstr = "tickle"  # string to use for arguments

        # waveform - parameter sweeps
        self.setattr_argument("att_tickle_db_list",
                              PYONValue([25]),
                              group=_argstr,
                              tooltip="Attenuation to be used for the urukul channel used for generating the tickle.")
        self.setattr_argument("ampl_tickle_pct_list",
                              PYONValue([50,]),
                              group=_argstr,
                              tooltip='Amplitude of tickle pulse.')

        self.setattr_argument("time_tickle_us_list",
                              PYONValue([1000]),
                              group=_argstr,
                              tooltip="Time for the total pulse (including pulse shape).")


        self.setattr_argument("type_pulse_shape",
                              EnumerationValue(list(available_pulse_shapes.keys()), default='square'),
                              group=_argstr,
                              tooltip="Pulse shape type to be used.")

        self.setattr_argument("freq_tickle_detuning_khz_list", Scannable(
            default=[
                ExplicitScan([0]),
                CenterScan(0., 10., 0.001, randomize=True),
                RangeScan(-10, 10, 26, randomize=True),
            ],
            global_min = -1000, global_max=1000, global_step=0.001,
            unit="kHz", scale=1, precision=6),
                              group=_argstr,
                              tooltip="Detuning from secular frequency of tickle pulse (in kHz) applied via the urukul dds.")

    def _build_arguments_intensity_servo(self):
        """
        Build arguements for intensity servo hold between shots
        """
        _argstr = 'intensity_servo_relock'
        self.setattr_argument('enable_servo_relock', BooleanValue(default=False), group=_argstr,
                              tooltip='Enables the servo to relock the intensity of the 729 beam after every shot')
        self.setattr_argument('time_servo_relock_us', NumberValue(default=2000, precision=3, step=1, min=1,
                                                                  max=10000, scale=1., unit='us'),
                              group=_argstr, tooltip='Length of time to let the servo relock before each shot')

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """

        ### Build Pulse Shaper
        self.dds_pulse_shaper = DDSPulseShaper(self, dds_target= self.dds_dipole.dds,
                                              ram_profile=self.profile_tickle_RAM,
                                              ram_addr_start=202, num_samples=250,
                                              ampl_max_pct=self.ampl_tickle_pct_list[0],
                                               pulse_shape=self.type_pulse_shape,
                                               phase_autoclear = 1)

        # run component preparation
        self._prepare_experiment_readout()

        freq_secular_ftw_list = self._prepare_experiment_ion_parameters()
        freq_tickle_detuning_ftw_list = self._prepare_experiment_tickle()

        # create experiment config
        self.config_experiment_list = create_experiment_config(

            # tickle sweeps
            freq_secular_ftw_list, freq_tickle_detuning_ftw_list,

            config_type=float, shuffle_config=True
        )

        # set timing for intensity servo between shots
        self.time_servo_relock_mu = self.core.seconds_to_mu(self.time_servo_relock_us * us)

    def _prepare_experiment_ion_parameters(self):
        """
        Prepare general ion parameters
        :return: list of carrier frequencies for cat and MS
        """
        self.freq_secular_ftw_list = [self.qubit.frequency_to_ftw(freq_secular_khz*kHz)
                                      for freq_secular_khz in self.freq_secular_khz_list]
        freq_secular_ftw_list = list(self.freq_secular_ftw_list)

        self.freq_carrier_ftw = self.qubit.frequency_to_ftw(self.freq_carrier_mhz * MHz)

        return freq_secular_ftw_list

    def _prepare_experiment_readout(self):
        """
        Prepare experiment values for state readout.
        """

        self.att_rap_mu = [att_to_mu(att_rap_db*dB) for att_rap_db in self.att_rap_db_list]
        self.rap_freq_dev_ftw_list = [self.qubit.frequency_to_ftw(freq_rap_dev_khz*kHz)
                                      for freq_rap_dev_khz in self.freq_rap_dev_khz_list]
        self.time_rap_mu_list = [self.core.seconds_to_mu(time_rap_us*us)
                                      for time_rap_us in self.time_rap_us_list]


        for idx, rap_profile in enumerate(self.freq_secular_khz_list):

            # attenuation register - readout (RAP): singlepasses set to default
            att_reg_readout_rap = 0x00000000 | (
                    (self.att_rap_mu[idx] << ((self.qubit.beam.chip_select - 4) * 8)) |
                    (self.qubit.att_singlepass0_default_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                    (self.qubit.att_singlepass1_default_mu << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                    (self.qubit.att_singlepass2_default_mu << ((self.qubit.singlepass2.chip_select - 4) * 8))
            )

            self.rap_subsequences.append(QubitRAP(
                self, ram_profile=idx+1, ram_addr_start=202, num_samples=250,
                ampl_max_pct=50, pulse_shape="blackman"
            ))

            self.att_reg_readout_rap_list.append(att_reg_readout_rap)

    def _prepare_experiment_tickle(self):
        """
        Prepare general experiment values for the tickle pulse.
        :return: tuple of (freq_tickle_detuning_hz_list)
        """
        # convert values to convenience units
        self.att_tickle_mu_list = [att_to_mu(att_tickle_db * dB) for att_tickle_db in self.att_tickle_db_list]
        freq_tickle_detuning_ftw_list = [self.dds_pulse_shaper.dds_target.frequency_to_ftw(freq_tickle_detuning_khz*kHz)
                                         for freq_tickle_detuning_khz in self.freq_tickle_detuning_khz_list]

        self.time_tickle_mu_list = [self.core.seconds_to_mu(time_heating_us * us) for
                                time_heating_us in self.time_tickle_us_list]

        return freq_tickle_detuning_ftw_list


    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                3)

    '''
    MAIN SEQUENCE
    '''

    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()
        self.dds_pulse_shaper.dds_target.sw.off()
        self.qubit.off()


    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # MAIN LOOP
        _loop_iter = 0
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:
                '''
                PREPARE & CONFIGURE
                '''
                # extract values from config list
                freq_secular_ftw = int32(config_vals[0])
                freq_tickle_detuning_ftw = int32(config_vals[1])

                '''
                BEGIN MAIN SEQUENCE
                '''
                self.core.break_realtime()  # add slack for execution
                delay_mu(125000)  # add even more slack lol
                time_tickle_mu = 0

                for idx in range(len(self.freq_secular_ftw_list)):
                    if freq_secular_ftw == self.freq_secular_ftw_list[idx]:

                        time_tickle_mu = self.time_tickle_mu_list[idx]

                        # configure tickle
                        self.dds_pulse_shaper.set_ampl_max_pct(self.ampl_tickle_pct_list[idx])
                        self.dds_pulse_shaper.sequence_initialize()
                        self.dds_pulse_shaper.dds_target.set_att_mu(self.att_tickle_mu_list[idx])
                        self.dds_pulse_shaper.dds_target.sw.off()
                        delay_mu(5000)
                        break

                '''
                Relock Intensity Servo
                '''
                if self.enable_servo_relock:
                    self.qubit.relock_intensity_servo(self.time_servo_relock_mu)

                self.ttl8.on()
                '''
                INITIALIZE ION STATE
                '''
                # initialize ion in S-1/2 state & SBC to ground state
                self.initialize_subsequence.run_dma()
                self.sidebandcool_subsequence.run_dma()
                delay_mu(8)

                # set tickle frequency/phases
                self.dds_pulse_shaper.dds_target.set_ftw(freq_secular_ftw + freq_tickle_detuning_ftw)
                self.dds_pulse_shaper.dds_target.set_pow(0)

                # set up config of shaped pulses to be fired for tickling
                # also sets up phase autoclear
                self.dds_pulse_shaper.configure_train(time_tickle_mu)
                self.dds_pulse_shaper.dds_target.cpld.io_update.pulse_mu(8)
                # for ururuk channel used for tickling keep RAM enabled but ensure we don't clear phase on io_update
                self.dds_pulse_shaper.dds_target.set_cfr1(ram_enable=1, phase_autoclear=0,
                                                          ram_destination=ad9910.RAM_DEST_ASF)
                self.dds_pulse_shaper.dds_target.cpld.io_update.pulse_mu(8)
                delay_mu(100)

                '''
                TICKLE PULSE
                '''
                self.dds_pulse_shaper.run_train_single()
                delay_mu(time_tickle_mu)

                '''
                READ OUT & STORE RESULTS
                '''
                self.pulse_readout_rap(freq_secular_ftw)

                # read out fluorescence & clean up loop
                self.readout_subsequence.run_dma()
                counts_res = self.readout_subsequence.fetch_count()

                # cleanup dds_pulse_shaper
                self.dds_pulse_shaper.sequence_cleanup()
                self.ttl8.off()

                # store results
                self.update_results(freq_secular_ftw,
                                    counts_res,
                                    freq_tickle_detuning_ftw,
                            )

                # check termination more frequently in case reps are low
                if _loop_iter % 100 == 0:
                    self.check_termination()
                _loop_iter += 1

            # rescue ion as needed & support graceful termination
            self.check_termination()

    @kernel(flags={"fast-math"})
    def pulse_readout_rap(self, sec_freq_ftw: TInt32) -> TNone:
        """
        Run a RAP readout pulse.
        """
        # set up relevant beam waveforms
        # find right sec freq
        rap_idx = 0
        for idx in range(len(self.freq_secular_ftw_list)):
            if sec_freq_ftw == self.freq_secular_ftw_list[idx]:
                rap_idx = idx
                break

        sec_freq_ftw = self.freq_secular_ftw_list[rap_idx]
        self.rap_subsequences[rap_idx].configure(self.time_rap_mu_list[rap_idx], self.freq_carrier_ftw - (sec_freq_ftw >> 1),
                                                 self.rap_freq_dev_ftw_list[rap_idx])
        delay_mu(50000)

        self.qubit.off()
        self.qubit.singlepass0_on()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()
        self.qubit.cpld.set_all_att_mu(self.att_reg_readout_rap_list[rap_idx])
        # run RAP readout pulse
        # run RAP turns on qubit
        self.rap_subsequences[rap_idx].run_rap(self.time_rap_mu_list[rap_idx])
