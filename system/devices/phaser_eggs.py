from numpy import int32, int64
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice
# todo: prepare - sample time, frame time
# todo: on/off func
# todo: att func
# todo: add conversion functions since artiq doesn't have them


class PhaserEGGS(LAXDevice):
    """
    Device: Phaser (EGGS)

    Use the Phaser AWG to create the EGGS RF
    """
    name = "phaser_eggs"
    core_device = ('phaser', 'phaser0')
    devices = {}
    kernel_invariants = {
        "t_sample_mu",
        "t_frame_mu",
        "t_output_delay_mu",
        "ftw_per_hz",
        "channel",
        "freq_center_hz",
        "freq_center_ftw",
        "phase_inherent_ch1_turns",
        "time_latency_ch1_system_ns"
    }

    def build_device(self):
        # set phaser sample/frame timings
        self.t_sample_mu =                  int64(40)
        self.t_frame_mu =                   int64(320)
        self.t_output_delay_mu =            int64(1953)
        # todo: max phaser sample rate

        # conversion factors
        # todo: fix - this is wrong-ish since everyone has a different ftw to hz conversion
        self.ftw_per_hz =                   (1 << 32) / 1e9

    def prepare_device(self):
        # alias both phaser output channels
        self.channel =                      self.phaser.channel

        # get frequency parameters
        self.freq_center_hz =               self.get_parameter('freq_center_mhz', group='eggs', override=False) * MHz
        self.freq_center_ftw =              self.frequency_to_ftw(self.freq_center_hz)

        # get phase delay parameters
        self.phase_inherent_ch1_turns =     self.get_parameter('phas_ch1_inherent_turns', group='eggs', override=False)
        self.time_latency_ch1_system_ns =   self.get_parameter('time_latency_ch1_system_ns', group='eggs', override=False)

    @kernel(flags={'fast-math'})
    def initialize_device(self):
        """
        Clear the DUC phase accumulators and sync the DAC.
        """
        # clear channel DUC phase accumulators
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_duc_cfg(clr_once=1)
        delay_mu(self.t_frame_mu)
        self.phaser.channel[1].set_duc_cfg(clr_once=1)
        delay_mu(self.t_frame_mu)

        # strobe update register for both DUCs
        self.phaser.duc_stb()
        delay_mu(self.t_frame_mu)

        # sync DAC for both channels
        self.phaser.dac_sync()
        # todo: set carrier frequency via DAC NCO frequency for both channels


    '''
    DUC Methods
    '''
    @kernel(flags={"fast-math"})
    def reset_duc_phase(self):
        """
        Disable amplitude and phase accumulator for all oscillators.
        """
        # synchronize to frame
        at_mu(self.phaser.get_next_frame_mu())

        # clear DUC phase accumulator for both channels
        self.phaser.channel[0].set_duc_cfg(clr_once=1)
        delay_mu(self.t_frame_mu)
        self.phaser.channel[1].set_duc_cfg(clr_once=1)
        delay_mu(self.t_frame_mu)

        # strobe update register for both DUCs
        self.phaser.duc_stb()


    '''
    Oscillator Methods
    '''

    @kernel(flags={"fast-math"})
    def disable_oscillators(self):
        """
        Set amplitude to 0 and keep phase accumulator cleared for all oscillators.
        # todo: document
        # note: this is different from reset_oscillators since it doesn't reset frequency, and persistently clears the phase accumulator
        """
        # synchronize to frame
        at_mu(self.phaser.get_next_frame_mu())

        # clear oscillator amplitudes
        for i in range(5):

            # do it for both channels
            with parallel:
                self.phaser.channel[0].oscillator[i].set_amplitude_phase(amplitude=0., clr=1)
                self.phaser.channel[1].oscillator[i].set_amplitude_phase(amplitude=0., clr=1)
                delay_mu(self.t_sample_mu)

    @kernel(flags={"fast-math"})
    def reset_oscillators(self):
        """
        Reset frequency and amplitude for all oscillators on both channels of the phaser.
        """
        for i in range(5):
            # synchronize to frame
            at_mu(self.phaser.get_next_frame_mu())

            # clear oscillator frequencies
            with parallel:
                self.phaser.channel[0].oscillator[i].set_frequency(0.)
                self.phaser.channel[1].oscillator[i].set_frequency(0.)
                delay_mu(self.t_sample_mu)

            # clear oscillator amplitudes
            with parallel:
                self.phaser.channel[0].oscillator[i].set_amplitude_phase(amplitude=0.)
                self.phaser.channel[1].oscillator[i].set_amplitude_phase(amplitude=0.)
                delay_mu(self.t_sample_mu)


    '''
    Helper Methods
    '''
    @portable(flags={"fast-math"})
    def amplitude_to_asf(self, amplitude: TFloat) -> TInt32:
        """
        todo: document
        """
        code = int32(round(amplitude * 0x7FFF))
        if code < 0 or code > 0x7FFF:
            raise ValueError("Error: Invalid fractional amplitude")
        return code

    @portable(flags={"fast-math"})
    def asf_to_amplitude(self, asf: TInt32) -> TFloat:
        """
        todo: document
        """
        return asf / float(0x7FFF)

    @portable(flags={"fast-math"})
    def turns_to_pow(self, turns: TFloat) -> TInt32:
        """
        todo: document
        """
        return int32(round(turns * 0x10000)) & int32(0xFFFF)

    @portable(flags={"fast-math"})
    def pow_to_turns(self, pow_: TInt32) -> TFloat:
        """
        todo: document
        """
        return pow_ / 0x10000

    @portable(flags={"fast-math"})
    def frequency_to_ftw(self, frequency: TFloat) -> TInt32:
        """
         todo: document
         """
        # todo: fix - this is wrong-ish since everyone has a different ftw to hz conversion
        return int32(round(self.ftw_per_hz * frequency))

    @portable(flags={"fast-math"})
    def ftw_to_frequency(self, ftw: TInt32) -> TFloat:
        """
         todo: document
         """
        # todo: fix - this is wrong-ish since everyone has a different ftw to hz conversion
        return ftw / self.ftw_per_hz
