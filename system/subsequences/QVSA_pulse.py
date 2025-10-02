from artiq.experiment import *
from numpy import array, zeros

from LAX_exp.language import *

from LAX_exp.system.objects.PulseShaper import available_pulse_shapes
from LAX_exp.system.objects.SpinEchoWizard import SpinEchoWizard
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper

# todo: migrate to SpinEchoWizardRDX
# todo: modernize


class QVSAPulse(LAXSubsequence):
    """
    Subsequence: QVSA Pulse

    Generate motional excitation using a multi-tone motional Raman transition (QVSA).
    Can be used to generate either a coherent or squeezed state.
    """
    name = 'qvsa_pulse'
    kernel_invariants = {
        # values
        'freq_qvsa_carrier_hz', 'freq_qvsa_secular_hz', 'att_qvsa_mu', 'waveform_qvsa_pulseshape_vals',

        # objects
        'spinecho_wizard', 'pulse_shaper'
    }

    def build_subsequence(self):
        # QVSA configuration - pulse
        self.setattr_argument("freq_qvsa_carrier_mhz",      NumberValue(default=80., precision=6, step=10., min=0., max=1000.), group='QVSA')
        self.setattr_argument("freq_qvsa_secular_khz",      NumberValue(default=1281., precision=3, step=1., min=0., max=5000.), group='QVSA')
        self.setattr_argument("time_qvsa_us",               NumberValue(default=1000, precision=2, step=500, min=0.04, max=100000000), group='QVSA')
        self.setattr_argument("ampl_qvsa_pct_config",       PYONValue([1., 1., 1.]), group='QVSA', tooltip="[rsb_qvsa_pct, bsb_qvsa_pct, carrier_qvsa_pct]")
        self.setattr_argument("phase_qvsa_turns_config",    PYONValue([0., 0., 0.]), group='QVSA', tooltip="[rsb_qvsa_turns, bsb_qvsa_turns, carrier_qvsa_turns]")
        self.setattr_argument("att_qvsa_db",                NumberValue(default=31.5, precision=1, step=0.5, min=0, max=31.5), group='QVSA')

        # QVSA configuration - pulse shaping
        self.setattr_argument("enable_pulse_shaping",           BooleanValue(default=False), group='QVSA.pulse_shaping')
        self.setattr_argument("type_pulse_shape",               EnumerationValue(list(available_pulse_shapes.keys()), default='sine_squared'), group='QVSA.pulse_shaping')
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=100, precision=1, step=100, min=0.2, max=100000), group='QVSA.pulse_shaping')
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=1000, precision=0, step=100, min=100, max=5000), group='QVSA.pulse_shaping')

        # get relevant devices
        self.setattr_device("core")
        self.setattr_device('phaser_eggs')

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizard(self)
        # set correct phase delays for field geometries (0.5 for osc_2 for dipole)
        self.pulse_shaper = PhaserPulseShaper(self, array([0., 0., 0.5, 0., 0.]))

    def prepare_subsequence(self):
        """
        Prepare values for speedy evaluation.
        """
        # validate inputs
        self._prepare_argument_checks()

        # convert QVSA values
        self.freq_qvsa_carrier_hz = self.freq_qvsa_carrier_mhz * MHz
        self.freq_qvsa_secular_hz = self.freq_qvsa_secular_khz * kHz
        # convert attenuation from dB to machine units
        self.att_qvsa_mu = att_to_mu(self.att_qvsa_db * dB)

        # configure waveform via pulse shaper & spin echo wizard
        self._prepare_waveform()

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # ensure phaser oscillator amplitudes are configured correctly
        if (type(self.ampl_qvsa_pct_config) is not list) or (len(self.ampl_qvsa_pct_config) != 3):
            raise ValueError("Invalid QVSA oscillator amplitude configuration."
                             "Must be of list [rsb_ampl_pct, bsb_ampl_pct, carrier_ampl_pct].")
        elif not all(0. <= val <= 100. for val in self.ampl_qvsa_pct_config):
            raise ValueError("Invalid QVSA oscillator amplitude. Must be in range [0., 100.].")
        elif sum(self.ampl_qvsa_pct_config) >= 100.:
            raise ValueError("Invalid QVSA oscillator amplitudes. Total must sum to <= 100.")

        # ensure phaser oscillator phases are configured correctly
        if (type(self.phase_qvsa_turns_config) is not list) or (len(self.phase_qvsa_turns_config) != 3):
            raise ValueError("Invalid QVSA oscillator phase configuration."
                             "Must be of list [rsb_phas_turns, bsb_phas_turns, carrier_phas_turns].")

        # ensure phaser output frequency falls within valid DUC bandwidth
        if abs(self.freq_qvsa_carrier_mhz * MHz - self.phaser_eggs.freq_center_hz) >= 200. * MHz:
            raise ValueError("Error: output frequencies outside +/- 300 MHz phaser DUC bandwidth.")

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the QVSA pulse.
        Uses SpinEchoWizard and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE WAVEFORM COMPILATION VIA SPINECHOWIZARD'''
        # create holding structures for EGGS pulse waveforms
        self._idx_waveform = 0

        # set up the spin echo wizard generally
        # note: time_pulse_us is amount of time for each block
        self.spinecho_wizard.time_pulse_us =                self.time_qvsa_us
        self.spinecho_wizard.enable_pulse_shaping =         self.enable_pulse_shaping
        self.spinecho_wizard.pulse_shape_blocks =           False
        self.spinecho_wizard.type_pulse_shape =             self.type_pulse_shape
        self.spinecho_wizard.time_pulse_shape_rolloff_us =  self.time_pulse_shape_rolloff_us
        self.spinecho_wizard.freq_pulse_shape_sample_khz =  self.freq_pulse_shape_sample_khz
        self.spinecho_wizard.enable_delay_spinecho =        False
        self.spinecho_wizard.time_delay_spinecho_us =       250

        '''DESIGN WAVEFORM SEQUENCE'''
        # create bare waveform block sequence
        _sequence_blocks = zeros((1, 3, 2), dtype=float)

        # set oscillator amplitudes
        _sequence_blocks[:, 0, 0] = self.ampl_qvsa_pct_config[0]
        _sequence_blocks[:, 1, 0] = self.ampl_qvsa_pct_config[1]
        _sequence_blocks[:, 2, 0] = self.ampl_qvsa_pct_config[2]

        # set oscillator phases (accounting for oscillator delay time)
        phase_bsb_update_delay_turns = (self.freq_qvsa_secular_khz * kHz) * (self.phaser_eggs.t_sample_mu * ns)
        _sequence_blocks[:, 0, 1] = self.phase_qvsa_turns_config[0]
        _sequence_blocks[:, 1, 1] = self.phase_qvsa_turns_config[1] + phase_bsb_update_delay_turns
        _sequence_blocks[:, 2, 1] = self.phase_qvsa_turns_config[2]

        # create waveform
        self.spinecho_wizard.sequence_blocks = _sequence_blocks
        self.spinecho_wizard.calculate_pulseshape()
        self.spinecho_wizard.compile_waveform()

        # get waveform data and store in holding structure
        self.waveform_qvsa_pulseshape_vals = self.spinecho_wizard.get_waveform()


    """
    KERNEL FUNCTIONS
    """
    @kernel(flags={"fast-math"})
    def initialize_subsequence(self) -> TNone:
        """
        Prepare the subsequence immediately before run.
        """
        ### PHASER INITIALIZATION ###
        # record phaser oscillator waveform
        # note: normally this would be encapsulated in phaser_record
        # e.g. in EGGSHeatingRDX-type exps
        delay_mu(1000000)  # add slack for recording DMA sequences (1000 us)
        self._idx_waveform = self.pulse_shaper.waveform_record(
            self.waveform_qvsa_pulseshape_vals[0],
            self.waveform_qvsa_pulseshape_vals[1],
            self.waveform_qvsa_pulseshape_vals[2]
        )
        self.core.break_realtime()

        # set maximum attenuations for phaser outputs to prevent leakage
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att_mu(0x00)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att_mu(0x00)
        delay_mu(25000)

        # configure phaser core hardware
        self.phaser_eggs.frequency_configure(self.freq_qvsa_carrier_hz,
                                             [-self.freq_qvsa_secular_hz, self.freq_qvsa_secular_hz, 0., 0., 0.],
                                             self.phaser_eggs.phase_inherent_ch1_turns)

    @kernel(flags={"fast-math"})
    def run_pulse(self) -> TInt64:
        """
        Run the main EGGS pulse together with supporting functionality.
        Returns:
            the start time of the phaser oscillator waveform.
            Useful to synchronize device operation.
        """
        # EGGS - START/SETUP
        self.phaser_eggs.phaser_setup(self.att_qvsa_mu, self.att_qvsa_mu)

        # EGGS - RUN
        # reset DUC phase to start DUC deterministically
        self.phaser_eggs.reset_duc_phase()
        # synchronize to next frame
        t_start_mu = self.phaser_eggs.get_next_frame_mu()
        at_mu(t_start_mu)
        self.pulse_shaper.waveform_playback(self._idx_waveform)

        # EGGS - STOP
        # stop all output & clean up hardware (e.g. eggs amp switches, RF integrator hold)
        # note: DOES unset attenuators (beware turn-on glitch if no filters/switches)
        self.phaser_eggs.phaser_stop()

        # return phaser osc start time (in case others want to sync)
        return t_start_mu

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        pass

