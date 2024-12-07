from numpy import int32, int64
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class PhaserEGGS(LAXDevice):
    """
    Device: Phaser (EGGS)

    Use the Phaser AWG to create the EGGS RF
    """
    name = "phaser_eggs"
    core_device = ('phaser', 'phaser1')
    devices = {
        'ch0_amp_sw':   'ttl8',     # switch AFTER EGGS CH0 amp
        'ch1_amp_sw':   'ttl9',     # switch AFTER EGGS CH1 amp
        'int_hold':     'ttl10'     # activate integrator hold on RF servo
    }
    kernel_invariants = {
        "t_sample_mu", "t_frame_mu", "t_output_delay_mu",
        "channel",
        "freq_center_hz", "phase_inherent_ch1_turns", "time_latency_ch1_system_ns",
        "time_rf_servo_holdoff_mu"
    }

    def build_device(self):
        # set phaser sample/frame timings
        self.t_sample_mu =          int64(40)
        self.t_frame_mu =           int64(320)
        self.t_output_delay_mu =    int64(1953)
        # todo: max phaser sample rate

    def prepare_device(self):
        # alias both phaser output channels
        self.channel = self.phaser.channel

        # get frequency parameters
        self.freq_center_hz = self.get_parameter('freq_center_mhz', group='eggs', override=False) * MHz

        # get phase delay parameters
        self.phase_inherent_ch1_turns =     self.get_parameter('phas_ch1_inherent_turns', group='eggs', override=False)
        self.time_latency_ch1_system_ns =   self.get_parameter('time_latency_ch1_system_ns', group='eggs', override=False)

        # get holdoff delay after phaser pulses to allow RF servo to re-lock
        self.time_rf_servo_holdoff_mu = self.get_parameter("time_rf_servo_holdoff_us", group="eggs",
                                                           conversion_function=us_to_mu)

    @kernel(flags={'fast-math'})
    def initialize_device(self) -> TNone:
        """
        Clear the DUC phase accumulators and sync the DAC.
        """
        # ensure INT HOLD on RF servo is deactivated
        self.int_hold.off()
        # ensure EGGS amp switches are closed
        self.ch0_amp_sw.off()
        self.ch1_amp_sw.off()

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

    @kernel(flags={'fast-math'})
    def cleanup_device(self) -> TNone:
        """
        Stop any residual output from the phaser.
        """
        self.core.break_realtime()
        self.phaser_stop()
        self.core.break_realtime()
        delay_mu(1000000)


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
        self.phaser.channel[0].set_att_mu(0x00)
        delay_mu(self.t_sample_mu)
        self.phaser.channel[1].set_att_mu(0x00)


    '''
    Setup/Prepare/Cleanup Methods
    '''
    @kernel(flags={"fast-math"})
    def phaser_setup(self, att_mu: TInt32) -> TNone:
        """
        Set up hardware to in preparation for an output pulse.
        Arguments:
            :param att_mu: attenuator value in machine units. 0x00 is 31.5 dB, 0xFF is 0 dB.
        """
        # EGGS - START/SETUP
        # set phaser attenuators - warning: creates turn on glitch
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_att_mu(att_mu)
        delay_mu(self.t_sample_mu)
        self.phaser.channel[1].set_att_mu(att_mu)

        # activate integrator hold
        self.int_hold.on()
        # add delay time after integrator hold to reduce effect of turn-on glitches
        delay_mu(self.time_rf_servo_holdoff_mu)

        # open phaser amp switches (add 2us delay for switches to fully open to prevent damage)
        with parallel:
            self.ch0_amp_sw.on()
            self.ch1_amp_sw.on()
        delay_mu(2000)

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

        # stop phaser amp switches - add extra delay b/c phaser pipeline pipeline latency
        delay_mu(5000)
        with parallel:
            self.ch0_amp_sw.off()
            self.ch1_amp_sw.off()

        # switch off EGGS attenuators to prevent phaser leakage
        delay_mu(self.t_sample_mu)
        self.phaser.channel[0].set_att_mu(0x00)
        delay_mu(self.t_sample_mu)
        self.phaser.channel[1].set_att_mu(0x00)

        # deactivate integrator hold - add delay time after EGGS pulse to allow RF servo to re-lock
        self.int_hold.off()
        delay_mu(self.time_rf_servo_holdoff_mu)
