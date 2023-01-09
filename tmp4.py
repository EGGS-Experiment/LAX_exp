from artiq.experiment import *
from artiq.coredevice.rtio import *
from artiq.coredevice.core import rtio_get_counter

#from numpy import int64


class Test3(EnvExperiment):
    """
    Test3
    Test the phaser0.
    """


    def build(self):
        self.setattr_device("core")
        self.setattr_device("phaser0")

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
        osc_num_list = [0, 1, 2, 3]
        ch_num = 0
        ch = self.phaser0.channel[ch_num]
        ph = self.phaser0
        self.core.break_realtime()


        # trf
        # ch.set_att(0 * dB)
        # ch.en_trf_out(rf=1, lo=0)
        # self.core.break_realtime()

        # duc
        ch.set_duc_frequency(10 * MHz)
        #ch.set_duc_cfg()
        ph.duc_stb()
        scb()


        '''
        *************ACTUAL*******************
        '''

        # maybe: have to write lo div sel (register 6, <<26)
        # pwd_out_buff
        # pwd_lo_div
        # pwd_tx_div
        # osc_freq_list = [-1.5, 1.5, 0.0, 0.0, 0.0]
        # osc_ampl_list = [0.49, 0.49, 0.0, 0.0, 0.0]
        #
        # for i in range(5):
        #     ch.oscillator[i].set_frequency(osc_freq_list[i] * MHz)
        #     ch.oscillator[i].set_amplitude_phase(amplitude=osc_ampl_list[i])
        #     self.core.break_realtime()
