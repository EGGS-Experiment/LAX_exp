from artiq.experiment import *
from artiq.coredevice import ad9910

from numpy import array, int32, int64, zeros
import numpy as np

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, Readout,
    RescueIon, QubitRAP
)

from LAX_exp.system.objects.PulseShaper import available_pulse_shapes
from LAX_exp.system.objects.dds_pulse_shaper import DDSPulseShaper

class MotionalCoherence(LAXExperiment, Experiment):
    """
    Experiment: Motional Coherence

    Perform two tickles seperated by a delay to test coherence of motional modes
    """
    name = 'Motional Coherence'
    kernel_invariants = {
        # subsequences & objects
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence',
        'rescue_subsequence', 'rap_subsequence', 'dds_pulse_shaper_tickle',

        # ion parameters
        'freq_secular_ftw',

        # hardware values - ramsey
        'enable_ramsey_delay',

        # hardware values - tickle
        'enable_tickle_pulse1', 'enable_tickle_pulse2', 'att_tickle_mu', 'time_tickle_mu',

        # hardware values - readout
        'freq_rap_center_ftw', 'freq_rap_dev_ftw', 'time_rap_mu',


        # configs
        'profile_729_SBC',
        'profile_729_RAP', 'profile_tickle_RAM',
        'att_reg_readout_rap',
        'config_experiment_list',

    }

    def build_experiment(self):

        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=50, precision=0, step=1, min=1, max=100000))

        # allocate relevant beam profiles
        self.profile_729_RAP = 0
        self.profile_729_SBC = 1

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
        self.setattr_device('ttl9')

        # set build arguments
        self._build_arguments_ion_parameters()
        self._build_arguments_ramsey()
        self._build_arguments_tickle_general()
        self._build_arguments_tickle1()
        self._build_arguments_tickle2()
        self._build_arguments_readout()

        # # instantiate RAP here since it relies on experiment arguments
        self.rap_subsequence = QubitRAP(
            self, ram_profile=self.profile_729_RAP, ram_addr_start=202, num_samples=250,
            ampl_max_pct=self.ampl_rap_pct, pulse_shape="blackman"
        )

    def _build_arguments_ion_parameters(self):
        """
        Build arguments for ion frequencies
        """
        _argstr = "ion_parameters"

        self.setattr_argument("freq_secular_khz", NumberValue(
            default=710,
            min=500, max=3000, step=0.001,
            unit="kHz", scale=1, precision=6),
                              group=_argstr,
                              tooltip="Secular frequency (in kHz) of the ion")

        self.setattr_argument("freq_carrier_mhz", NumberValue(
            default=100.41,
            min=60.,
            max=400,
            step=1,
            unit="MHz",
            scale=1,
            precision=6
        ),
          group=_argstr,
          tooltip="Carrier frequency of the ion.\n"
                  "Note: this is applied via the main doublepass DDS.\n")


    def _build_arguments_readout(self):
        """
        Build arguments for readout pulse.
        """
        # RAP-based readout
        self.setattr_argument("att_rap_db",
                              NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, unit="dB", scale=1.),
                              group="read.RAP")
        self.setattr_argument("ampl_rap_pct",
                              NumberValue(default=50., precision=3, step=5, min=1, max=50, unit="%", scale=1.),
                              group="read.RAP")
        self.setattr_argument("freq_rap_dev_khz",
                              NumberValue(default=72., precision=2, step=0.01, min=1, max=1e4, unit="kHz", scale=1.),
                              group='read.RAP')
        self.setattr_argument("time_rap_us",
                              NumberValue(default=400., precision=3, min=1, max=1e7, step=1, unit="us", scale=1.),
                              group="read.RAP")

    def _build_arguments_tickle_general(self):
        """
        Build core sweep arguments for the tickle pulses.
        """
        _argstr = "tickle_general"  # string to use for arguments

        # waveform - parameter sweeps
        self.setattr_argument("att_tickle_db",
                              NumberValue(default=31., precision=1, step=0.5, min=0., max=31.5, unit="dB", scale=1.),
                              group=_argstr,
                              tooltip="Attenuation to be used for the urukul channel used for generating the tickle.")
        self.setattr_argument("ampl_tickle_pct",
                              NumberValue(default=50., precision=2, min=0., max=50., unit="%", scale=1.),
                              group=_argstr,
                              tooltip='Amplitude of tickle pulse.')

        self.setattr_argument("time_heating_us",
                              NumberValue(default=50, precision=2, step=500, min=0.04, max=10000000, unit="us",
                                          scale=1.),
                              group=_argstr,
                              tooltip="Time for the total pulse (including pulse shape).")

        # waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping", BooleanValue(default=False),
                              group=_argstr,
                              tooltip="Applies pulse shaping to the edges of the tickle pulse.")
        self.setattr_argument("type_pulse_shape",
                              EnumerationValue(list(available_pulse_shapes.keys()), default='sine_squared'),
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

    def _build_arguments_tickle1(self):

        _argstr = "tickle1"
        self.setattr_argument('enable_tickle_pulse1',
                              BooleanValue(default=True),
                              tooltip='turn on first tickle pulse',
                              group= _argstr)

    def _build_arguments_tickle2(self):

        _argstr = "tickle2"
        self.setattr_argument('enable_tickle_pulse2',
                              BooleanValue(default=True),
                              tooltip='turn on first tickle pulse',
                             group=_argstr)


        self.setattr_argument("phase_tickle2_offset_turns_list", Scannable(
            default=[
                ExplicitScan([0.]),
                RangeScan(0, 1.0, 26, randomize=True),
            ],
            global_min=0.0, global_max=1.0, global_step=1,
            unit="turns", scale=1, precision=5),
                              group=_argstr,
                              tooltip="Phase of second tickle pulse (in turns) "
                                      "relative to the first tickle pulse"
                                      "applied via the urukul dds.")

    def _build_arguments_ramsey(self):
        _argstr = 'ramsey'
        # ramsey delay between tickl1 and tickle2
        self.setattr_argument("enable_ramsey_delay", BooleanValue(default=False), group=_argstr,
                              tooltip="Enables a Ramsey delay between the 1st and 2nd tickle pulses. "
                                      "Useful for doing motional coherence tests.")
        self.setattr_argument("time_ramsey_delay_us_list", Scannable(
            default=[
                ExplicitScan([100]),
                RangeScan(0, 500, 50, randomize=True), ],
            global_min=1, global_max=100000, global_step=1,
            unit="us", scale=1, precision=5),
                              group=_argstr,
                              tooltip="Ramsey delay time between 1st and 2nd tickle pulses.")

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        ### Build Pulse Shaper
        if self.enable_pulse_shaping:
            pulse_shape = self.type_pulse_shape
        else:
            pulse_shape = 'square'
        self.dds_pulse_shaper_tickle = DDSPulseShaper(self, dds_target= self.dds_dipole.dds,
                                              ram_profile=self.profile_tickle_RAM,
                                              ram_addr_start=202, num_samples=200,
                                              ampl_max_pct=self.ampl_tickle_pct,
                                               pulse_shape=pulse_shape,
                                               phase_autoclear = 1,
                                            external_switch=self.ttl9)

        # run component preparation
        self._prepare_experiment_ion_parameters()
        self._prepare_experiment_readout()
        freq_tickle_detuning_ftw_list = self._prepare_experiment_tickle_general()
        phase_tickle2_offset_pow_list = self._prepare_experiment_tickle2()
        time_ramsey_delay_mu_list = self._prepare_experiment_ramsey()

        # create experiment config
        self.config_experiment_list = create_experiment_config(
            time_ramsey_delay_mu_list,
            # tickle sweeps
            freq_tickle_detuning_ftw_list, phase_tickle2_offset_pow_list,
            config_type=float, shuffle_config=True
        )

    def _prepare_experiment_ion_parameters(self):
        """
        Prepare general ion parameters
        :return: list of carrier frequencies for cat and MS
        """
        self.freq_secular_ftw = self.qubit.frequency_to_ftw(self.freq_secular_khz*kHz)

        self.freq_carrier_ftw = self.qubit.frequency_to_ftw(self.freq_carrier_mhz * MHz)

    def _prepare_experiment_ramsey(self):
        # inter-cat ramsey delay
        if self.enable_ramsey_delay:
            time_ramsey_delay_mu_list = [self.core.seconds_to_mu(time_delay_us * us)
                                         for time_delay_us in self.time_ramsey_delay_us_list]
        else:
            time_ramsey_delay_mu_list = [0]

        return time_ramsey_delay_mu_list

    def _prepare_experiment_readout(self):
        """
        Prepare experiment values for state readout.
        """
        # prepare RAP arguments
        self.freq_rap_center_ftw = self.freq_carrier_ftw - (self.freq_secular_ftw >> 1)
        self.freq_rap_dev_ftw = self.qubit.frequency_to_ftw(self.freq_rap_dev_khz * kHz)
        self.time_rap_mu = self.core.seconds_to_mu(self.time_rap_us * us)


        # attenuation register - readout (RAP): singlepasses set to default
        self.att_reg_readout_rap = 0x00000000 | (
                (att_to_mu(self.att_rap_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass0_default_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass1_default_mu << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass2_default_mu << ((self.qubit.singlepass2.chip_select - 4) * 8)))

    def _prepare_experiment_tickle_general(self):
        """
        Prepare general experiment values for the tickle pulse.
        :return: tuple of (freq_tickle_detuning_hz_list, phase_tickle_list)
        """
        # convert values to convenience units
        self.att_tickle_mu = att_to_mu(self.att_tickle_db * dB)
        freq_tickle_detuning_ftw_list = [self.dds_pulse_shaper_tickle.dds_targets[0].frequency_to_ftw(freq_tickle_detuning_khz*kHz)
                                         for freq_tickle_detuning_khz in self.freq_tickle_detuning_khz_list]
        self.time_tickle_mu = self.core.seconds_to_mu(self.time_heating_us * us)

        return freq_tickle_detuning_ftw_list

    def _prepare_experiment_tickle2(self):
        phase_tickle2_offset_pow_list = [self.dds_pulse_shaper_tickle.dds_targets[0].turns_to_pow(phase_tickle_turns) for
                                 phase_tickle_turns in self.phase_tickle2_offset_turns_list]
        return phase_tickle2_offset_pow_list

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                4)

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

        # configure RAP pulse
        self.rap_subsequence.configure(self.time_rap_mu,
                                       self.freq_rap_center_ftw,
                                       self.freq_rap_dev_ftw)
        delay_mu(50000)

        # configure tickle
        self.dds_pulse_shaper_tickle.sequence_initialize()
        self.dds_pulse_shaper_tickle.dds_targets[0].set_att_mu(self.att_tickle_mu)
        self.dds_pulse_shaper_tickle.dds_targets[0].sw.off()
        delay_mu(50000)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # MAIN LOOP
        _loop_iter = 0
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:
                '''
                PREPARE & CONFIGURE
                '''
                time_ramsey_delay_mu = int64(config_vals[0])
                freq_tickle_detuning_ftw = int32(config_vals[1])
                phase_tickle2_offset_pow = int32(config_vals[2])

                '''
                BEGIN MAIN SEQUENCE
                '''
                self.core.break_realtime()  # add slack for execution
                delay_mu(125000)  # add even more slack lol

                '''
                INITIALIZE ION STATE
                '''
                # initialize ion in S-1/2 state & SBC to ground state
                self.initialize_subsequence.run_dma()
                self.sidebandcool_subsequence.run_dma()

                # set tickle frequency/phases
                self.dds_pulse_shaper_tickle.dds_targets[0].set_ftw(
                    self.freq_secular_ftw + freq_tickle_detuning_ftw)
                self.dds_pulse_shaper_tickle.dds_targets[0].set_pow(0)
                # set up config of shaped pulses to be fired for tickling, also sets up phase autoclear
                time_actual_tickle_list_mu = self.dds_pulse_shaper_tickle.configure_train_all_dds(
                    [self.time_tickle_mu])
                time_actual_tickle_mu = time_actual_tickle_list_mu[0]

                # # set tickle frequency/phases
                ref_time_mu = self.clear_tickle_phase_accumulators(
                    self.dds_pulse_shaper_tickle
                )

                '''
                TICKLE PULSE 1 
                '''
                if self.enable_tickle_pulse1:
                    time_tickle1_start_mu = self.dds_pulse_shaper_tickle.run_train_all_dds()
                '''
                RAMSEY DELAY
                '''
                if self.enable_ramsey_delay:
                    delay_mu(time_ramsey_delay_mu)
                '''
                TICKLE PULSE 2
                '''
                if self.enable_tickle_pulse2:
                    self.dds_pulse_shaper_tickle.dds_targets[0].set_pow(phase_tickle2_offset_pow)
                    time_tickle2_start_mu = self.dds_pulse_shaper_tickle.run_train_all_dds()

                '''
                READ OUT & STORE RESULTS
                '''

                self.pulse_readout_rap()
                # read out fluorescence & clean up loop
                self.readout_subsequence.run_dma()
                counts_res = self.readout_subsequence.fetch_count()
                # cleanup dds_pulse_shaper_tickle
                self.dds_pulse_shaper_tickle.sequence_cleanup()

                # store results
                self.update_results(freq_tickle_detuning_ftw,
                                    counts_res,
                                    phase_tickle2_offset_pow,
                                    time_ramsey_delay_mu,
                                    )

                # check termination more frequently in case reps are low
                if _loop_iter % 100 == 0:
                    self.check_termination()
                _loop_iter += 1

            # rescue ion as needed & support graceful termination
            self.check_termination()


    @kernel(flags={"fast-math"})
    def clear_tickle_phase_accumulators(self,
                                            dds_pulse_shaper_tickle) -> TInt64:
        # set cfr1 so we clear phases of all urukul0 channels on next io_update
        # tickle dds has already had phase autoclear flag set high in drg setup

        ref_time_mu = (now_mu() + 8) & ~7
        at_mu(ref_time_mu)
        dds_pulse_shaper_tickle.dds_targets[0].cpld.io_update.pulse_mu(8)
        # for ururuk channel used for tickling keep RAM enabled but ensure we don't clear phase on io_update
        dds_pulse_shaper_tickle.dds_targets[0].set_cfr1(ram_enable=1, phase_autoclear=0,
                                                             ram_destination=ad9910.RAM_DEST_ASF)
        dds_pulse_shaper_tickle.dds_targets[0].cpld.io_update.pulse_mu(8)

        return ref_time_mu

    @kernel(flags={"fast-math"})
    def pulse_readout_rap(self) -> TNone:
        """
        Run a RAP readout pulse.
        """
        # set up relevant beam waveforms
        self.qubit.off()
        self.qubit.singlepass0_on()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()
        self.qubit.singlepass0.set_mu(self.qubit.freq_singlepass0_default_ftw,
                                      asf=self.qubit.ampl_singlepass0_default_asf, pow_=0,
                                      profile=self.profile_729_RAP,
                                      phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
        delay_mu(2000)
        self.qubit.cpld.set_all_att_mu(self.att_reg_readout_rap)
        # run RAP readout pulse
        # run RAP turns on qubit
        self.rap_subsequence.run_rap(self.time_rap_mu)

    def analyze_experiment(self):
        pass