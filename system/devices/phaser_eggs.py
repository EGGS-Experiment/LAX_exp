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
        "t_sample_mu", "t_frame_mu", "t_output_delay_mu",
        "ftw_per_hz", "channel",
        "freq_center_hz", "phase_inherent_ch1_turns", "time_latency_ch1_system_ns"
    }

    def build_device(self):
        # set phaser sample/frame timings
        self.t_sample_mu =          int64(40)
        self.t_frame_mu =           int64(320)
        self.t_output_delay_mu =    int64(1953)
        # todo: max phaser sample rate

        # conversion factors
        # todo: fix - this is wrong-ish since everyone has a different ftw to hz conversion
        self.ftw_per_hz =           (1 << 32) / 1e9

    def prepare_device(self):
        # alias both phaser output channels
        self.channel =                      self.phaser.channel

        # get frequency parameters
        self.freq_center_hz =               self.get_parameter('freq_center_mhz', group='eggs', override=False) * MHz

        # get phase delay parameters
        self.phase_inherent_ch1_turns =     self.get_parameter('phas_ch1_inherent_turns', group='eggs', override=False)
        self.time_latency_ch1_system_ns =   self.get_parameter('time_latency_ch1_system_ns', group='eggs', override=False)

    @kernel(flags={'fast-math'})
    def initialize_device(self) -> TNone:
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

    @kernel(flags={'fast-math'})
    def cleanup_device(self) -> TNone:
        """
        Stop any residual output from the phaser.
        """
        self.core.break_realtime()
        delay_mu(1000000)
        print("debug: cleanup phaser")
        self.core.break_realtime()

        self.reset_oscillators()


    '''
    DUC Methods
    '''
    @kernel(flags={"fast-math"})
    def reset_duc_phase(self) -> TNone:
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
    def reset_oscillators(self) -> TNone:
        """
        Reset frequency, amplitude, and phase accumulators for all oscillators on both channels of the phaser.
        Maximum attenuation is set to minimize leakage of downstream phaser components.
        This function synchronizes to the frame for each oscillator reset.
        """
        # clear oscillators
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
                self.phaser.channel[0].oscillator[i].set_amplitude_phase(amplitude=0., clr=1)
                self.phaser.channel[1].oscillator[i].set_amplitude_phase(amplitude=0., clr=1)
                delay_mu(self.t_sample_mu)

        # set max attenuations for phaser outputs to reduce effect of internal noise
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_att(31.5 * dB)
        delay_mu(self.t_sample_mu)
        self.phaser.channel[1].set_att(31.5 * dB)

    @kernel(flags={"fast-math"})
    def phaser_stop(self) -> TNone:
        """
        Stop the phaser quickly.
        Set maximum attenuation to prevent output leakage.
        """
        # disable eggs phaser output
        with parallel:
            self.phaser.channel[0].oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            self.phaser.channel[1].oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser.channel[0].oscillator[1].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            self.phaser.channel[1].oscillator[1].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser.channel[0].oscillator[2].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            self.phaser.channel[1].oscillator[2].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser.channel[0].oscillator[3].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            self.phaser.channel[1].oscillator[3].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser.channel[0].oscillator[4].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            self.phaser.channel[1].oscillator[4].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            delay_mu(self.t_sample_mu)

        # switch off EGGS attenuators to prevent leakage
        # delay_mu(self.t_sample_mu)
        # self.phaser.channel[0].set_att(31.5 * dB)
        # delay_mu(self.t_sample_mu)
        # self.phaser.channel[1].set_att(31.5 * dB)
