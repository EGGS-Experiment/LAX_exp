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
        # get devices
        ch0 = self.phaser0.channel[0]
        ch1 = self.phaser0.channel[1]

        # initialize
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()
        self.phaser0.init(debug=False)
        delay(1*ms)

        # set attenuations
        ch0.set_att(0*dB)
        delay(1*ms)
        ch1.set_att(0*dB)
        delay(1*ms)

        # set DUC
        ch0.set_duc_frequency(0*MHz)
        delay(1*ms)
        ch1.set_duc_frequency(0*MHz)
        delay(1*ms)
        ch0.set_duc_phase(.0)
        delay(1*ms)
        ch1.set_duc_phase(.0)
        delay(1*ms)
        ch0.set_duc_cfg(clr_once=1)
        delay(1*ms)
        ch1.set_duc_cfg(clr_once=1)
        delay(.1*ms)
        self.phaser0.duc_stb()
        delay(1*ms)

        t_stamp = self.get_frame_timestamp()
        at_mu(t_stamp + 100000)

        ch0.oscillator[0].set_frequency(100 * MHz)
        delay_mu(40 * 2)
        ch0.oscillator[0].set_amplitude_phase(amplitude=0.18, phase=0.)
        delay(10*ms)

        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

        # # loop
        # for _ in range(10000):
        #
        #     # ensure updates are aligned to a frame
        #     t_stamp = self.get_frame_timestamp()
        #     at_mu(t_stamp + 100000)
        #
        #     # output a waveform
        #     ch0.oscillator[0].set_frequency(150 * MHz)
        #     delay_mu(40*2)
        #     ch0.oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=1)
        #     delay_mu(40*10)
        #     ch0.oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=0)
        #     delay_mu(40*2)
        #     ch0.oscillator[0].set_amplitude_phase(amplitude=.8)
        #     delay_mu(40*25)
        #     ch0.oscillator[0].set_amplitude_phase(amplitude=.8, phase=0.5)
        #     delay_mu(40*25)
        #     ch0.oscillator[0].set_amplitude_phase(amplitude=.0)
        #
        #     # add break
        #     self.core.wait_until_mu(now_mu())
        #     # print("make jitter")
        #     self.core.break_realtime()
        #     delay(10*ms)
