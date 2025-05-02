from numpy import int64
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class PhaserEGGS(LAXDevice):
    """
    Device: Phaser (EGGS)

    Use the Phaser AWG to create the EGGS RF.
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
        "time_rf_servo_holdoff_mu", "time_switch_holdoff_mu"
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

        # get holdoff delay after setting TTL switches
        self.time_switch_holdoff_mu = self.get_parameter("time_switch_holdoff_us", group="eggs",
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
    Setup/Prepare/Cleanup Methods
    '''
    @kernel(flags={"fast-math"})
    def frequency_configure(self, carrier_freq_hz: TFloat, osc_freq_hz_list: TList(TFloat),
                            phase_ch1_offset_turns: TFloat) -> TNone:
        """
        Configure oscillator frequencies on phaser for EGGS.
        Arguments:
            carrier_freq_hz: the output center frequency (in Hz).
                center frequency is set identically on CH0 and CH1.
            osc_freq_hz_list: list of frequencies to set each oscillator.
                oscillators are set identically on CH0 and CH1.
            phase_ch1_offset_turns: the phase offset for CH1 relative to CH0 (in turns).
        """
        '''
        CALCULATE PHASE DELAYS
        '''
        # calculate phase delays between CH0 and CH1, accounting for the relative CH1 latency
        phase_ch1_turns = phase_ch1_offset_turns + (carrier_freq_hz * self.time_latency_ch1_system_ns * ns)
        delay_mu(100000) # add slack

        '''
        SET CARRIER FREQUENCY
        '''
        # set carrier offset frequency via the DUC
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_duc_frequency(carrier_freq_hz - self.freq_center_hz)
        delay_mu(self.t_frame_mu)
        self.phaser.channel[1].set_duc_frequency(carrier_freq_hz - self.freq_center_hz)
        delay_mu(self.t_frame_mu)
        # strobe updates for both channels
        self.phaser.duc_stb()

        # set DUC phase delay to compensate for inter-channel latency
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[1].set_duc_phase(phase_ch1_turns)
        self.phaser.duc_stb()

        '''
        SET OSCILLATOR (i.e. sideband) FREQUENCIES
        '''
        # synchronize to frame, then set oscillator frequencies
        at_mu(self.phaser.get_next_frame_mu())
        with parallel:
            self.phaser.channel[0].oscillator[0].set_frequency(osc_freq_hz_list[0])
            self.phaser.channel[1].oscillator[0].set_frequency(osc_freq_hz_list[0])
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser.channel[0].oscillator[1].set_frequency(osc_freq_hz_list[1])
            self.phaser.channel[1].oscillator[1].set_frequency(osc_freq_hz_list[1])
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser.channel[0].oscillator[2].set_frequency(osc_freq_hz_list[2])
            self.phaser.channel[1].oscillator[2].set_frequency(osc_freq_hz_list[2])
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser.channel[0].oscillator[3].set_frequency(osc_freq_hz_list[3])
            self.phaser.channel[1].oscillator[3].set_frequency(osc_freq_hz_list[3])
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser.channel[0].oscillator[4].set_frequency(osc_freq_hz_list[4])
            self.phaser.channel[1].oscillator[4].set_frequency(osc_freq_hz_list[4])
            delay_mu(self.t_sample_mu)

    @kernel(flags={"fast-math"})
    def phaser_setup(self, att_mu_ch0: TInt32 = 0, att_mu_ch1: TInt32 = 0) -> TNone:
        """
        Set up hardware in preparation for an output pulse.
        Arguments:
            :param att_mu_ch0: phaser CH0 attenuator value in machine units. 0x00 is 31.5 dB, 0xFF is 0 dB.
            :param att_mu_ch1: phaser CH1 attenuator value in machine units. 0x00 is 31.5 dB, 0xFF is 0 dB.
        """
        # EGGS - START/SETUP
        # set phaser attenuators - warning: creates turn on glitch
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_att_mu(att_mu_ch0)
        delay_mu(self.t_sample_mu)
        self.phaser.channel[1].set_att_mu(att_mu_ch1)

        with parallel:
            # open phaser amp switches (add 25us delay for switches to fully open to prevent damage)
            self.ch0_amp_sw.on()
            self.ch1_amp_sw.on()
            # activate integrator hold
            self.int_hold.on()

        # add delay time after integrator hold to reduce effect of turn-on glitches
        delay_mu(self.time_rf_servo_holdoff_mu)
        # note: need 25us delay b/c max rise/fall time of ZSW2-272VDHR+ switches
        delay_mu(self.time_switch_holdoff_mu)

    @kernel(flags={"fast-math"})
    def phaser_stop(self) -> TNone:
        """
        Stop the phaser quickly.
        Set maximum attenuation to prevent output leakage.
        Can be used following an EGGS pulse to stop the phaser.
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
        delay_mu(2500)
        with parallel:
            self.ch0_amp_sw.off()
            self.ch1_amp_sw.off()
            # deactivate integrator hold
            self.int_hold.off()

        # add delay time after EGGS pulse to allow RF servo to re-lock
        delay_mu(self.time_rf_servo_holdoff_mu)

        # switch off EGGS attenuators to prevent phaser leakage
        delay_mu(self.t_sample_mu)
        self.phaser.channel[0].set_att_mu(0x00)
        delay_mu(self.t_sample_mu)
        self.phaser.channel[1].set_att_mu(0x00)


