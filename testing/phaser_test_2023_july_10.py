
from artiq.experiment import *
from artiq.coredevice.rtio import *
from artiq.coredevice.core import rtio_get_counter

from numpy import int64


class PhaserTest(EnvExperiment):
    """
    Phaser Test
    Test the phaser.
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("phaser0")

        # todo: select channel
        # todo: select att
        # todo: select on/off
        # todo: select NCO freq
        # todo: select DUC freq
        # todo: select carrier
        # todo: select oscillator amplitudes & frequency

    def prepare(self):
        # alias break_realtime() since i'm lazy
        self.scb = self.core.break_realtime

        # phaser timings
        self.time_sample_mu =                                           int64(40)
        self.time_frame_mu =                                            self.phaser0.t_frame

    @kernel
    def get_frame_timestamp(self):
        """Performs a Phaser register read and records the exact frame timing in self.frame_tstamp.
        For deterministic pulses, changes must be phase-aligned with the phaser frames
        """
        rtio_output(self.phaser0.channel_base << 8, 0)  # read some phaser reg (board id here)
        delay_mu(int64((8*4*10)))
        frame_tstamp = rtio_input_timestamp(rtio_get_counter()+0xffffff, self.phaser0.channel_base)
        delay(10*ms)
        return frame_tstamp

    @kernel
    def run(self):
        # setup
        self.core.reset()

        # scb = self.core.break_realtime
        # ch_num = 0
        # ch = self.phaser0.channel[ch_num]
        # ch1 = self.phaser0.channel[1]
        # self.core.break_realtime()
        #
        # '''
        # *************INIT*******************
        # '''
        # self.phaser0.init(debug=True)
        # scb()
        #
        # # nco stuff
        # at_mu(self.phaser0.get_next_frame_mu())
        # ch.set_nco_frequency((-217.083495 - 0) * MHz)
        # at_mu(self.phaser0.get_next_frame_mu())
        # ch1.set_nco_frequency((-217.083495 - 0) * MHz)
        # self.core.break_realtime()
        #
        # at_mu(self.phaser0.get_next_frame_mu())
        # ch.set_nco_phase(0.)
        # delay_mu(self.time_sample_mu)
        # ch1.set_nco_phase(0.)
        # self.core.break_realtime()
        #
        # at_mu(self.phaser0.get_next_frame_mu())
        # self.phaser0.dac_sync()
        # self.core.break_realtime()
        #
        #
        # # attenuators
        # at_mu(self.phaser0.get_next_frame_mu())
        # self.phaser0.channel[0].set_att(10 * dB)
        # at_mu(self.phaser0.get_next_frame_mu())
        # self.phaser0.channel[1].set_att(10 * dB)
        #
        # # trf output
        # at_mu(self.phaser0.get_next_frame_mu())
        # ch.en_trf_out(rf=1, lo=0)
        # at_mu(self.phaser0.get_next_frame_mu())
        # ch1.en_trf_out(rf=1, lo=0)
        # self.core.break_realtime()
        #
        #
        # # duc
        # at_mu(self.phaser0.get_next_frame_mu())
        # ch.set_duc_frequency(0. * MHz)
        # at_mu(self.phaser0.get_next_frame_mu())
        # ch1.set_duc_frequency(0. * MHz)
        # self.core.break_realtime()
        #
        # at_mu(self.phaser0.get_next_frame_mu())
        # ch.set_duc_cfg()
        # at_mu(self.phaser0.get_next_frame_mu())
        # ch1.set_duc_cfg()
        # self.core.break_realtime()
        #
        # at_mu(self.phaser0.get_next_frame_mu())
        # self.phaser0.duc_stb()
        # self.core.break_realtime()
        #
        #
        # '''
        # *************DMA*******************
        # '''
        #
        # # prepare DMA
        # self._prepare_dma()
        # self.core.break_realtime()
        #
        # # retrieve DMA handle
        # handle = self.core_dma.get_handle('_PHASER_SET')
        # self.core.break_realtime()
        #
        #
        # '''
        # *************ACTUAL*******************
        # '''
        #
        # # align to frame and clear DUC phase and synchronize DAC
        # at_mu(self.phaser0.get_next_frame_mu())
        # self.phaser0.dac_sync()
        #
        # at_mu(self.phaser0.get_next_frame_mu())
        # self.phaser0.channel[0].set_duc_cfg(clr=1)
        # self.phaser0.channel[1].set_duc_cfg(clr=1)
        # at_mu(self.phaser0.get_next_frame_mu())
        # self.phaser0.duc_stb()
        #
        # # measure frame timestamp
        # self.phaser0.measure_frame_timestamp()
        #
        #
        # # run DMA
        # # self.core_dma.playback_handle(handle)
        # # self.core.break_realtime()
        #
        # # tmp remove
        #
        # at_mu(self.phaser0.get_next_frame_mu())
        # self.phaser0.channel[0].set_duc_cfg(clr=0)
        # self.phaser0.channel[1].set_duc_cfg(clr=0)
        # at_mu(self.phaser0.get_next_frame_mu())
        # self.phaser0.duc_stb()

        # clear each oscillator amplitude and phase
        at_mu(self.phaser0.get_next_frame_mu())
        # with parallel:
        self.phaser0.channel[0].oscillator[0].set_amplitude_phase(amplitude=0.999, phase=0.0, clr=1)
        self.phaser0.channel[1].oscillator[0].set_amplitude_phase(amplitude=0.999, phase=0.0, clr=1)
        # delay_mu(self.time_sample_mu)

        # set each oscillator amplitude and phase
        at_mu(self.phaser0.get_next_frame_mu())
        # with parallel:
        self.phaser0.channel[0].oscillator[0].set_amplitude_phase(amplitude=0.999, phase=0.0, clr=0)
        self.phaser0.channel[1].oscillator[0].set_amplitude_phase(amplitude=0.999, phase=0.585, clr=0)
        delay_mu(self.time_sample_mu)

        # tmp remove

        # maybe: have to write lo div sel (register 6, <<26)
        # pwd_out_buff
        # pwd_lo_div
        # pwd_tx_div

        # osc_freq_list = [1.696, -1.696, 0.0, 0.0, 0.0]
        # osc_ampl_list = [0.499, 0.499, 0.0, 0.0, 0.0]
        #
        # for i in range(5):
        #     at_mu(self.phaser0.get_next_frame_mu())
        #     ch.oscillator[i].set_frequency(osc_freq_list[i] * MHz)
        #     at_mu(self.phaser0.get_next_frame_mu())
        #     ch.oscillator[i].set_amplitude_phase(amplitude=osc_ampl_list[i], clr=0)
        #     self.core.break_realtime()

        # # kill all output
        # at_mu(self.phaser0.get_next_frame_mu())
        # self.phaser0.channel[0].en_trf_out(rf=1, lo=0)
        # self.phaser0.channel[1].en_trf_out(rf=1, lo=0)
        # self.core.break_realtime()


    @kernel(flags={"fast-math"})
    def _prepare_dma(self):
        self.core.break_realtime()

        # alias channels for conveniece
        ch0 = self.phaser0.channel[0]
        ch1 = self.phaser0.channel[1]

        # set up values
        osc_freq_list = [0.0, 0.0, 0.0, 0.0, 0.0]
        osc_ampl_list = [0.999, 0.0, 0.0, 0.0, 0.0]

        osc1_phase_list = [0.0, 0.0, 0.0, 0.0, 0.0]
        osc2_phase_list = [0.0, 0.0, 0.0, 0.0, 0.0]
        self.core.break_realtime()
        delay_mu(100000)
        self.core.break_realtime()

        # set up frequencies
        for i in range(5):
            at_mu(self.phaser0.get_next_frame_mu())
            with parallel:
                ch0.oscillator[i].set_frequency(osc_freq_list[i] * MHz)
                ch1.oscillator[i].set_frequency(osc_freq_list[i] * MHz)

            at_mu(self.phaser0.get_next_frame_mu())
            with parallel:
                ch0.oscillator[i].set_amplitude_phase(amplitude=osc_ampl_list[i], phase=osc1_phase_list[0], clr=1)
                ch1.oscillator[i].set_amplitude_phase(amplitude=osc_ampl_list[i], phase=osc2_phase_list[0], clr=1)

        # set nco phase and dac sync
        #         at_mu(self.phaser0.get_next_frame_mu())
        #         ch.set_nco_phase(0.)
        #         delay_mu(self.time_sample_mu)
        #         ch1.set_nco_phase(0.)
        #         self.core.break_realtime()
        #
        #         at_mu(self.phaser0.get_next_frame_mu())
        #         self.phaser0.dac_sync()
        #         self.core.break_realtime()

        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.dac_sync()

        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.measure_frame_timestamp()

        # record dma sequence
        with self.core_dma.record('_PHASER_SET'):

            # # align to frame and clear DUC phase
            # at_mu(self.phaser0.get_next_frame_mu())
            # ch0.set_duc_cfg(clr=1)
            # ch1.set_duc_cfg(clr=1)
            # at_mu(self.phaser0.get_next_frame_mu())
            # self.phaser0.duc_stb()
            #
            # # clear oscillator phase accumulator
            # at_mu(self.phaser0.get_next_frame_mu())
            # # with parallel:
            # ch0.oscillator[0].set_amplitude_phase(amplitude=osc_ampl_list[0], phase=osc1_phase_list[0], clr=1)
            # ch1.oscillator[0].set_amplitude_phase(amplitude=osc_ampl_list[0], phase=osc2_phase_list[0], clr=1)

            # align to frame and clear DUC phase
            at_mu(self.phaser0.get_next_frame_mu())
            ch0.set_duc_cfg(clr=0)
            ch1.set_duc_cfg(clr=0)
            at_mu(self.phaser0.get_next_frame_mu())
            self.phaser0.duc_stb()

            # set each oscillator amplitude and phase
            at_mu(self.phaser0.get_next_frame_mu())
            with parallel:
                ch0.oscillator[0].set_amplitude_phase(amplitude=osc_ampl_list[0], phase=osc1_phase_list[0], clr=0)
                ch1.oscillator[0].set_amplitude_phase(amplitude=osc_ampl_list[0], phase=osc2_phase_list[0], clr=0)
            delay_mu(self.time_sample_mu)

            # # set osc 1
            # with parallel:
            #     ch0.oscillator[1].set_amplitude_phase(amplitude=osc_ampl_list[1], phase=osc1_phase_list[1], clr=0)
            #     ch1.oscillator[1].set_amplitude_phase(amplitude=osc_ampl_list[1], phase=osc2_phase_list[1], clr=0)
            # delay_mu(self.time_sample_mu)
            #
            # # set osc 2
            # ch0.oscillator[2].set_amplitude_phase(amplitude=osc_ampl_list[2], phase=osc1_phase_list[2], clr=0)
            # delay_mu(self.time_sample_mu)
            # ch1.oscillator[2].set_amplitude_phase(amplitude=osc_ampl_list[2], phase=osc2_phase_list[2], clr=0)
            # delay_mu(self.time_sample_mu)
            #
            # # set osc 2
            # ch0.oscillator[3].set_amplitude_phase(amplitude=osc_ampl_list[3], phase=osc1_phase_list[3], clr=0)
            # delay_mu(self.time_sample_mu)
            # ch1.oscillator[3].set_amplitude_phase(amplitude=osc_ampl_list[3], phase=osc2_phase_list[3], clr=0)
            # delay_mu(self.time_sample_mu)
            #
            # # set osc 4
            # ch0.oscillator[4].set_amplitude_phase(amplitude=osc_ampl_list[4], phase=osc1_phase_list[4], clr=0)
            # delay_mu(self.time_sample_mu)
            # ch1.oscillator[4].set_amplitude_phase(amplitude=osc_ampl_list[4], phase=osc2_phase_list[4], clr=0)
            # delay_mu(self.time_sample_mu)
