
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
        self.setattr_device("phaser0")

        # todo: select channel
        # todo: select att
        # todo: select on/off
        # todo: select NCO freq
        # todo: select DUC freq
        # todo: select carrier
        # todo: select oscillator amplitudes & frequency

    def prepare(self):
        self.scb = self.core.break_realtime

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

        scb = self.core.break_realtime
        ch_num = 0
        ch = self.phaser0.channel[ch_num]
        self.core.break_realtime()

        '''
        *************INIT*******************
        '''
        self.phaser0.init(debug=True)
        scb()

        # nco stuff
        at_mu(self.phaser0.get_next_frame_mu())
        ch.set_nco_frequency((-217.083495 - 0) * MHz)
        # ch.set_nco_frequency(-302.083495* MHz)
        self.core.break_realtime()

        at_mu(self.phaser0.get_next_frame_mu())
        ch.set_nco_phase(0.)
        self.core.break_realtime()

        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.dac_sync()
        self.core.break_realtime()


        # trf
        at_mu(self.phaser0.get_next_frame_mu())
        ch.set_att(3 * dB)
        self.core.break_realtime()

        at_mu(self.phaser0.get_next_frame_mu())
        ch.en_trf_out(rf=1, lo=0)
        self.core.break_realtime()


        # duc
        at_mu(self.phaser0.get_next_frame_mu())
        ch.set_duc_frequency(-2 * MHz)
        self.core.break_realtime()

        at_mu(self.phaser0.get_next_frame_mu())
        ch.set_duc_cfg()
        self.core.break_realtime()

        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.duc_stb()
        self.core.break_realtime()


        '''
        *************ACTUAL*******************
        '''

        # maybe: have to write lo div sel (register 6, <<26)
        # pwd_out_buff
        # pwd_lo_div
        # pwd_tx_div
        osc_freq_list = [1.693, -1.693, 0.0, 0.0, 0.0]
        osc_ampl_list = [0.499, 0.499, 0.0, 0.0, 0.0]

        for i in range(5):
            at_mu(self.phaser0.get_next_frame_mu())
            ch.oscillator[i].set_frequency(osc_freq_list[i] * MHz)
            at_mu(self.phaser0.get_next_frame_mu())
            ch.oscillator[i].set_amplitude_phase(amplitude=osc_ampl_list[i], clr=0)
            self.core.break_realtime()

        # # kill all output
        # at_mu(self.phaser0.get_next_frame_mu())
        # ch.en_trf_out(rf=1, lo=0)
        # self.core.break_realtime()
