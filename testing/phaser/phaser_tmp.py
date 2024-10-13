import numpy as np
from artiq.experiment import *


class PhaserSet(EnvExperiment):
    """
    Utility: Phaser Set

    Set the phaser output for fast user testing.
    """
    name = 'Phaser Set'

    def build(self):
        # get scheduler
        self.setattr_device("scheduler")

        # global arguments
        self.setattr_argument("repetitions",                    NumberValue(default=5, precision=0, step=1, min=1, max=10000))
        self.setattr_argument("cleanup",                        BooleanValue(default=False), group='global')
        self.setattr_argument("freq_carrier_mhz",               NumberValue(default=82, precision=5, step=1, min=1, max=10000), group='global')
        self.setattr_argument("time_pulse_ms",                  NumberValue(default=10000., precision=5, step=1, min=0.000001, max=100000), group='global')
        self.setattr_argument("time_reset_ms",                  NumberValue(default=2., precision=5, step=1, min=0.000001, max=100000), group='global')
        # self.setattr_argument("clear_dac_phase_accumulator", BooleanValue(default=False), group='global')
        self.setattr_argument("synchronize_oscillator_updates", BooleanValue(default=True), group='global')

        # channel-specific arguments
        self.setattr_argument("att_ch0_db",                     NumberValue(default=0., precision=1, step=0.5, min=0, max=31.5), group='channels')
        self.setattr_argument("att_ch1_db",                     NumberValue(default=0., precision=1, step=0.5, min=0, max=31.5), group='channels')
        self.setattr_argument("time_latency_ch1_system_ns",     NumberValue(default=0.17, precision=3, step=0.1, min=-1.0, max=1.0), group='channels')
        self.setattr_argument("apply_ch1_latency_on_oscillators", BooleanValue(default=True), group='channels')
        self.setattr_argument("phase_inherent_ch1_turns",       NumberValue(default=0.07, precision=3, step=0.1, min=-1.0, max=1.0), group='channels')

        # oscillators - channel 0
        # note: oscillator waveforms are in the format [freq_khz, ampl_pct, phas_turns]
        self.setattr_argument("waveform_ch0_osc0",  PYONValue(default=[-1000., 49., 0.]), group='oscillators_ch0')
        self.setattr_argument("waveform_ch0_osc1",  PYONValue(default=[1000., 49., 0.]), group='oscillators_ch0')
        self.setattr_argument("waveform_ch0_osc2",  PYONValue(default=[0., 0., 0.]), group='oscillators_ch0')
        self.setattr_argument("waveform_ch0_osc3",  PYONValue(default=[0., 0., 0.]), group='oscillators_ch0')
        self.setattr_argument("waveform_ch0_osc4",  PYONValue(default=[0., 0., 0.]), group='oscillators_ch0')

        # oscillators - channel 1
        # note: oscillator waveforms are in the format [freq_khz, ampl_pct, phas_turns]
        self.setattr_argument("waveform_ch1_osc0",  PYONValue(default=[-1000., 49., 0.]), group='oscillators_ch1')
        self.setattr_argument("waveform_ch1_osc1",  PYONValue(default=[1000., 49., 0.]), group='oscillators_ch1')
        self.setattr_argument("waveform_ch1_osc2",  PYONValue(default=[0., 0., 0.]), group='oscillators_ch1')
        self.setattr_argument("waveform_ch1_osc3",  PYONValue(default=[0., 0., 0.]), group='oscillators_ch1')
        self.setattr_argument("waveform_ch1_osc4",  PYONValue(default=[0., 0., 0.]), group='oscillators_ch1')

        # build constants
        self._build_constants()

    def _build_constants(self):
        """
        Instantiate fixed/constant objects/variables.
        """
        # get relevant devices
        self.setattr_device('core')
        self.setattr_device('core_dma')
        self.setattr_device('phaser0')
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')

        # hardware values
        self.t_sample_mu =          np.int64(40)
        self.t_frame_mu =           np.int64(320)
        self.time_output_delay_mu = np.int64(1953)

        # preallocate DAC34H84 0x1F register config to speed up
        # writes when clearing the phase
        self.dac_register_1f =      np.int32(0)

        # preallocate phase holders
        self.phase_ch1_turns =      np.float(0)

        self.phase_ch0_osc0 =       np.float(0)
        self.phase_ch0_osc1 =       np.float(0)
        self.phase_ch0_osc2 =       np.float(0)
        self.phase_ch0_osc3 =       np.float(0)
        self.phase_ch0_osc4 =       np.float(0)

        self.phase_ch1_osc0 =       np.float(0)
        self.phase_ch1_osc1 =       np.float(0)
        self.phase_ch1_osc2 =       np.float(0)
        self.phase_ch1_osc3 =       np.float(0)
        self.phase_ch1_osc4 =       np.float(0)


    def prepare(self):
        # set relevant values for phaser initialization
        self.time_phaser_sample_mu = np.int64(40)

    @kernel
    def run(self):
        self.core.reset()

        '''
        *************PHASER*******************
        '''
        # initialize phaser
        self.phaser0.init(debug=True)


        '''
        *************DAC*******************
        '''
        # set DAC NCO frequency to center output at 85 MHz exactly
        # note: TRF372017 freq is 302.083918 MHz => DAC NCO should be 217.083918 MHz
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_nco_frequency((-217.083495) * MHz)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_nco_frequency((-217.083495) * MHz)


        # clear DAC NCO phase
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_nco_phase(0.)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_nco_phase(0.)

        # sync DAC for both channels
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.dac_sync()


        '''
        *************DUC (DIGITAL UPCONVERTERS)*******************
        '''
        # set channel DUC frequencies to 0
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_frequency(0 * MHz)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_duc_frequency(0 * MHz)

        # clear channel DUC phase accumulators
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_cfg(clr_once=1)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_duc_cfg(clr_once=1)

        # strobe update register for both DUCs
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.duc_stb()


        '''
        *************OSCILLATORS*******************
        '''
        # reset oscillator frequency and amplitude, and keep phase accumulator persistently cleared
        # note: this has to happen before TRF or attenuator adjustment to ensure channel outputs are 0
        for i in range(5):
            # clear channel 0 oscillator
            at_mu(self.phaser0.get_next_frame_mu())
            self.phaser0.channel[0].oscillator[i].set_frequency(0 * MHz)
            delay_mu(self.time_phaser_sample_mu)
            self.phaser0.channel[0].oscillator[i].set_amplitude_phase(amplitude=0., clr=1)

            # clear channel 1 oscillator
            at_mu(self.phaser0.get_next_frame_mu())
            self.phaser0.channel[1].oscillator[i].set_frequency(0 * MHz)
            delay_mu(self.time_phaser_sample_mu)
            self.phaser0.channel[1].oscillator[i].set_amplitude_phase(amplitude=0., clr=1)


        '''
        *************TRF (UPCONVERTER)*******************
        '''
        # enable outputs for both channels
        # note: want to leave trf outputs persistently enabled since phase relation
        # between channels can change after adjusting the TRF
        # note: no need to set TRF frequency here since we already do this in device_db
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].en_trf_out(rf=1, lo=0)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].en_trf_out(rf=1, lo=0)


        '''
        *************ATTENUATORS*******************
        '''
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_att(31.5 * dB)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_att(31.5 * dB)

        self.core.break_realtime()




        # center @ 85
        at_mu(self.phaser0.get_next_frame_mu())
        # self.phaser0.channel[0].set_nco_frequency(-217.083495 * MHz)
        self.phaser0.channel[0].set_nco_frequency(0. * MHz)
        at_mu(self.phaser0.get_next_frame_mu())
        # self.phaser0.channel[1].set_nco_frequency(-217.083495 * MHz)
        self.phaser0.channel[1].set_nco_frequency(0. * MHz)
        # sync DAC for both channels
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.dac_sync()

        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_frequency(0. * MHz)
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.duc_stb()
        self.core.break_realtime()



        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_att(0. * dB)

        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].oscillator[0].set_frequency(1.0 * MHz)
        # at_mu(self.phaser0.get_next_frame_mu())
        # self.phaser0.channel[0].oscillator[1].set_frequency(2 * MHz)
        # at_mu(self.phaser0.get_next_frame_mu())
        # self.phaser0.channel[0].oscillator[2].set_frequency(0 * MHz)

        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].oscillator[0].set_amplitude_phase(amplitude=0.999, phase=0.0, clr=0)
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].oscillator[1].set_amplitude_phase(amplitude=0.0, phase=0.0, clr=0)
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].oscillator[2].set_amplitude_phase(amplitude=0.0, phase=0.0, clr=0)
        self.core.break_realtime()

