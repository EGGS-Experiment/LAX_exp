from numpy import int32, int64
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice
# todo: add function to dynamically change TRF frequency


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
        "t_sample_mu", "t_frame_mu", "t_output_delay_mu", "t_output_delay_mu", "t_output_delay_mu",
        "channel",
        "freq_center_hz", "phase_inherent_ch1_turns", "time_latency_ch1_system_ns",
        "time_phaser_holdoff_mu"
    }

    def build_device(self):
        # set phaser sample/frame timings
        self.t_sample_mu =          int64(40)   # oscillator update sample period
        self.t_frame_mu =           int64(320)  # fastlink frame period
        self.t_output_delay_mu =    int64(1953) # measured pipeline latency relative to a TTL
        # forced delay between setup/stop and sequence finish (6 frame periods; deliberately close to t_output_delay_mu)
        self.t_setup_delay_mu =     int64(1920)

    def prepare_device(self):
        # alias both phaser output channels
        self.channel = self.phaser.channel

        # get frequency parameters
        self.freq_center_hz = self.get_parameter('freq_center_mhz', group='devices.phaser', override=False) * MHz

        # get phase delay parameters
        self.phase_inherent_ch1_turns =     self.get_parameter('phas_ch1_inherent_turns', group='devices.phaser.ch1', override=False)
        self.time_latency_ch1_system_ns =   self.get_parameter('time_latency_ch1_system_ns', group='devices.phaser.ch1', override=False)

        # holdoff delay before & after phaser pulses (e.g. to allow RF servo to re-lock, to account for switch rise times)
        self.time_phaser_holdoff_mu = self.get_parameter("time_phaser_holdoff_us", group="devices.phaser", conversion_function=us_to_mu)
        # ensure delays are multiples of phaser frame period
        self.time_phaser_holdoff_mu = int64(round(self.time_phaser_holdoff_mu / self.t_frame_mu) * self.t_frame_mu)

        self.switch_delay_time_mu = int64(8)

    @kernel(flags={'fast-math'})
    def initialize_device(self) -> TNone:
        """
        Clear the DUC phase accumulators and sync the DAC.
        """
        self.amp_sws_off()

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
        self.phaser_stop()


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
        :param carrier_freq_hz: the output center frequency (in Hz).
            center frequency is set identically on CH0 and CH1.
        :param osc_freq_hz_list: list of frequencies to set each oscillator.
            oscillators are set identically on CH0 and CH1.
        :param phase_ch1_offset_turns: the phase offset for CH1 relative to CH0 (in turns).
        """
        '''
        CALCULATE PHASE DELAYS
        '''
        # calculate phase delays between CH0 and CH1, accounting for the relative CH1 latency
        phase_ch1_turns = phase_ch1_offset_turns + (carrier_freq_hz * self.time_latency_ch1_system_ns * ns)
        delay_mu(99840) # add slack (320ns multiple)

        '''
        SET CARRIER FREQUENCY
        '''
        # set carrier offset frequency via the DUC
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_duc_frequency(carrier_freq_hz - self.freq_center_hz)
        # delay_mu(self.t_frame_mu)
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[1].set_duc_frequency(carrier_freq_hz - self.freq_center_hz)
        # delay_mu(self.t_frame_mu)
        at_mu(self.phaser.get_next_frame_mu())
        # strobe updates for both channels
        self.phaser.duc_stb()

        # set DUC phase delay to compensate for inter-channel latency
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[1].set_duc_phase(phase_ch1_turns)
        self.phaser.duc_stb()

        # # get fiducial start time
        # time_start_mu = self.phaser.get_next_frame_mu()
        #
        # # set carrier offset frequency via the DUC and strobe updates for both channels
        # at_mu(time_start_mu)
        # self.phaser.channel[0].set_duc_frequency(carrier_freq_hz - self.freq_center_hz)
        # at_mu(time_start_mu + self.t_frame_mu)
        # self.phaser.channel[1].set_duc_frequency(carrier_freq_hz - self.freq_center_hz)
        # at_mu(time_start_mu + 2 * self.t_frame_mu)
        # self.phaser.duc_stb()
        #
        # # set DUC phase delay to compensate for inter-channel latency
        # at_mu(time_start_mu + 3 * self.t_frame_mu)
        # self.phaser.channel[1].set_duc_phase(phase_ch1_turns)
        # self.phaser.duc_stb()

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
        :param att_mu_ch0: phaser CH0 attenuator value in machine units. 0x00 is 31.5 dB, 0xFF is 0 dB.
        :param att_mu_ch1: phaser CH1 attenuator value in machine units. 0x00 is 31.5 dB, 0xFF is 0 dB.
        """
        # EGGS - START/SETUP
        # set phaser attenuators - warning: creates turn on glitch
        #   must do while switches are closed (ideally do first for maximum margin)
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_att_mu(att_mu_ch0)
        delay_mu(self.t_frame_mu)
        self.phaser.channel[1].set_att_mu(att_mu_ch1)

        # add extra delay to ensure turn-on glitches are suppressed by switches
        delay_mu(self.t_setup_delay_mu - 16)  # 6 frame periods (minus coarse RTIO between switches)

        self.amp_sws_on()

        # add holdoff delay to account for integrator hold, switch rise/fall times (25us for ZSW2-272VDHR+)
        delay_mu(self.time_phaser_holdoff_mu)

    @kernel(flags={"fast-math"})
    def phaser_stop(self) -> TNone:
        """
        Stop the phaser quickly.
        Set maximum attenuation to prevent output leakage.
        Can be used following an EGGS pulse to stop the phaser.
        """
        # disable phaser oscillator output
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

        # add delay for oscillator updates to account for pipeline latency
        delay_mu(self.t_setup_delay_mu - 16)  # 8 frame periods (minus coarse RTIO between switches)
        # stop phaser amp switches & deactivate integrator hold

        self.amp_sws_off()

        # add delay time after EGGS pulse to allow RF servo to re-lock
        delay_mu(self.time_phaser_holdoff_mu)

        # switch off EGGS attenuators to prevent phaser leakage
        # note: want these to happen after holdoff to ensure attenuator video feedthru
        #   stopped by switches (ideally do last)
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_att_mu(0x00)
        delay_mu(self.t_frame_mu)
        self.phaser.channel[1].set_att_mu(0x00)

    @kernel(flags={"fast-math"})
    def phase_osc_clear(self) -> TNone:
        """
        Clear oscillator phase accumulators.
        Note: this function does not account for the 40ns-ish phaser sample periods.
            Just sit and hope the phase delay isn't too bad lol.
        """
        # clear oscillator phases
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

        # reenable oscillator phase accumulators
        with parallel:
            self.phaser.channel[0].oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=0)
            self.phaser.channel[1].oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=0)
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser.channel[0].oscillator[1].set_amplitude_phase(amplitude=0., phase=0., clr=0)
            self.phaser.channel[1].oscillator[1].set_amplitude_phase(amplitude=0., phase=0., clr=0)
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser.channel[0].oscillator[2].set_amplitude_phase(amplitude=0., phase=0., clr=0)
            self.phaser.channel[1].oscillator[2].set_amplitude_phase(amplitude=0., phase=0., clr=0)
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser.channel[0].oscillator[3].set_amplitude_phase(amplitude=0., phase=0., clr=0)
            self.phaser.channel[1].oscillator[3].set_amplitude_phase(amplitude=0., phase=0., clr=0)
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser.channel[0].oscillator[4].set_amplitude_phase(amplitude=0., phase=0., clr=0)
            self.phaser.channel[1].oscillator[4].set_amplitude_phase(amplitude=0., phase=0., clr=0)
            delay_mu(self.t_sample_mu)

    @kernel(flags={'fast-math'})
    def ch0_amp_sw_on(self) -> TNone:
        self.ch0_amp_sw.on()
        delay_mu(self.switch_delay_time_mu)

    @kernel(flags={'fast-math'})
    def ch0_amp_sw_off(self) -> TNone:
        self.ch0_amp_sw.off()
        delay_mu(self.switch_delay_time_mu)

    @kernel(flags={'fast-math'})
    def ch1_amp_sw_on(self) -> TNone:
        self.ch1_amp_sw.on()
        delay_mu(self.switch_delay_time_mu)

    @kernel(flags={'fast-math'})
    def ch1_amp_sw_off(self) -> TNone:
        self.ch1_amp_sw.off()
        delay_mu(self.switch_delay_time_mu)

    @kernel(flags={'fast-math'})
    def int_hold_off(self) -> TNone:
        self.int_hold.off()
        delay_mu(self.switch_delay_time_mu)

    @kernel(flags={'fast-math'})
    def int_hold_on(self) -> TNone:
        self.int_hold.on()
        delay_mu(self.switch_delay_time_mu)

    @kernel(flags={'fast-math'})
    def amp_sws_on(self) -> TNone:
        self.ch0_amp_sw_on()
        self.ch1_amp_sw_on()
        self.int_hold_on()
        delay_mu(self.time_phaser_holdoff_mu)

    @kernel(flags={'fast-math'})
    def amp_sws_off(self) -> TNone:
        self.ch0_amp_sw_off()
        self.ch1_amp_sw_off()
        self.int_hold_off()
        delay_mu(self.time_phaser_holdoff_mu)


    '''
    Conversion functions
    '''
    @portable(flags={"fast-math"})
    def nco_frequency_to_ftw(self, freq_hz: TFloat = 0.) -> TInt32:
        """
        Return the 32-bit frequency tuning word corresponding to the given
        frequency for the DAC34H84 NCO.
        :param freq_hz: freq_hz: NCO frequency in Hz, passband from [-400, 400] MHz.
        :returns: 32-bit frequency tuning word (in Hz).
        """
        ftw = int32(round(freq_hz * ((1 << 30)/(250 * MHz))))
        if ftw < 0x0 or ftw > 0xFFFFFFFF:
            raise ValueError("Invalid NCO frequency, must be in range [-500, 500] MHz")
        return ftw

    @portable(flags={"fast-math"})
    def duc_frequency_to_ftw(self, freq_hz: TFloat = 0.) -> TInt32:
        """
        Return the 32-bit frequency tuning word corresponding to the given
        frequency for the Phaser's DUC.
        :param freq_hz: freq_hz: NCO frequency in Hz, passband from [-250, 250] MHz.
        :returns: 32-bit frequency tuning word (in Hz).
        """
        ftw = int32(round(freq_hz * ((1 << 30)/(125 * MHz))))
        if ftw < 0x0 or ftw > 0xFFFFFFFF:
            raise ValueError("Invalid DUC frequency, must be in range [-250, 250] MHz")
        return ftw

    @portable(flags={"fast-math"})
    def osc_frequency_to_ftw(self, freq_hz: TFloat = 0.) -> TInt32:
        """
        Return the 32-bit frequency tuning word corresponding to the given
        frequency for a Phaser MultiDDS oscillator.
        :param freq_hz: freq_hz: NCO frequency in Hz, passband from [-250, 250] MHz.
        :returns: 32-bit frequency tuning word (in Hz).
        """
        ftw = int32(round(freq_hz * ((1 << 30)/(6.25 * MHz))))
        if ftw < 0x0 or ftw > 0xFFFFFFFF:
            raise ValueError("Invalid oscillator frequency, must be in range [-12.5, 12.5] MHz")
        return ftw

