from artiq.experiment import *
import LAX_exp.experiments.ion_spectrum_analyzer.IonSpectrumAnalyzer as IonSpectrumAnalyzer


class IonSpectrumAnalyzerRamsey(IonSpectrumAnalyzer.IonSpectrumAnalyzer):
    """
    Experiment: Ion Spectrum Analyzer Ramsey

    Apply the Ion Spectrum Analyzer with some delay time to conduct a Ramsey-style experiment.
    """
    name = 'Ion Spectrum Analyzer Ramsey'


    def build_experiment(self):
        # run regular ISA build
        super().build_experiment()

        # ISA ramsey arguments
        self.setattr_argument("enable_ISA_antisqueezing_ramsey",            BooleanValue(default=False), group='ISA.antisqueezing')
        self.setattr_argument("time_ISA_antisqueezing_ramsey_us",           NumberValue(default=50., ndecimals=3, step=100, min=0, max=1000000), group='ISA.antisqueezing')

    def prepare_experiment(self):
        super().prepare_experiment()

        # configure ISA ramsey behavior
        self.time_ISA_antisqueezing_ramsey_us =     self.core.seconds_to_mu(self.time_ISA_antisqueezing_ramsey_us * us)
        if self.enable_ISA_antisqueezing_ramsey:    self.ISA_ramsey = self._phaser_run_ramsey
        else:                                       self.ISA_ramsey = self.noop()


    '''
    HARDWARE
    '''
    @kernel(flags={"fast-math"})
    def phaser_run_ISA_antisqueezing(self, ampl_rsb_frac: TFloat, ampl_bsb_frac: TFloat, ampl_dd_frac: TFloat):
        """
        Activate phaser channel outputs for EGGS heating.
        Sets the same RSB, BSB, and dynamical decoupling amplitudes for both channels.
        Adjusts the BSB phase by a set value after a set time (should be calculated in _prepare_squeezing).
        Arguments:
            ampl_rsb_frac   (float) : the red sideband amplitude (as a decimal fraction).
            ampl_bsb_frac   (float) : the blue sideband amplitude (as a decimal fraction).
            ampl_dd_frac    (float) : the dynamical decoupling amplitude (as a decimal fraction).
        """
        '''
        ISA - INTRINSIC SQUEEZING
        '''
        # set oscillator 0 (RSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_ch0_osc0, clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_ch1_osc0, clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_ch0_osc1, clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_ch1_osc1, clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_ch0_osc2, clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_ch1_osc2, clr=0)

        # heat for first half (before antisqueezing)
        delay_mu(self.time_ISA_antisqueeze_mu)

        # do customizable ramsey
        self.ISA_ramsey()

        '''
        ISA - INTRINSIC ANTISQUEEZING
        '''
        # note: no need to account for 40ns delay time between oscillator updates since the oscillator phases are NOT cleared
        # adjust oscillator 0 (BSB) phase for antisqueezing
        with parallel:
            self.ttl8.on()
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=self.ampl_ISA_antisqueezing_rsb_frac, phase=self.phase_ch0_osc0 + self.phase_ISA_antisqueezing_rsb_turns, clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=self.ampl_ISA_antisqueezing_rsb_frac, phase=self.phase_ch1_osc0 + self.phase_ISA_antisqueezing_rsb_turns, clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # adjust oscillator 1 (BSB) phase for antisqueezing
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=self.ampl_ISA_antisqueezing_bsb_frac, phase=self.phase_ch0_osc1 + self.phase_ISA_antisqueezing_bsb_turns, clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=self.ampl_ISA_antisqueezing_bsb_frac, phase=self.phase_ch1_osc1 + self.phase_ISA_antisqueezing_bsb_turns, clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # turn off oscillator 2 (carrier) during antisqueezing
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=self.ampl_ISA_antisqueezing_dipole_frac, phase=self.phase_ch0_osc2 + self.phase_ISA_antisqueezing_dipole_turns, clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=self.ampl_ISA_antisqueezing_dipole_frac, phase=self.phase_ch1_osc2 + self.phase_ISA_antisqueezing_dipole_turns, clr=0)

        # heat for second half
        delay_mu(self.time_ISA_antisqueeze_mu)
        self.ttl8.off()

    @kernel(flags={"fast-math"})
    def _phaser_run_ramsey(self):
        """
        Set up hardware for the reverse ISA ramsey pulse (i.e. ISA antisqueezing).
        Note: we DON'T clear the oscillators here to keep the reverse pulse phase-coherent with the forward pulse.
        """
        # CLEAR OSCILLATORS
        # set oscillator 0 (RSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=0., phase=self.phase_ch0_osc0, clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=0., phase=self.phase_ch1_osc0, clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=0., phase=self.phase_ch0_osc1, clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=0., phase=self.phase_ch1_osc1, clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=0., phase=self.phase_ch0_osc2, clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=0., phase=self.phase_ch1_osc2, clr=0)

        # delay ramsey time
        delay_mu(self.time_ISA_antisqueezing_ramsey_us)
