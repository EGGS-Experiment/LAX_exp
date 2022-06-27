from artiq.experiment import *
from artiq.coredevice.rtio import *
from artiq.coredevice.core import rtio_get_counter


class PhaserTest(EnvExperiment):
    """
    Phaser Test
    Test the phaser.
    """


    def build(self):
        self.setattr_device("core")
        self.setattr_device("phaser")

    @kernel
    def get_frame_timestamp(self):
        """Performs a Phaser register read and records the exact frame timing in self.frame_tstamp.
        For deterministic pulses, changes must be phase-aligned with the phaser frames
        """
        rtio_output(self.phaser.channel_base << 8, 0)  # read some phaser reg (board id here)
        delay_mu(int64((8*4*10)))
        frame_tstamp = rtio_input_timestamp(rtio_get_counter()+0xffffff, self.phaser.channel_base)
        delay(10*ms)
        return frame_tstamp

    @kernel
    def run(self):
        ch0 = self.phaser.channel[0]
        ch1 = self.phaser.channel[1]
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()
        self.phaser.init(debug=False)
        delay(1*ms)
        ch0.set_att(0*dB)
        delay(1*ms)
        ch1.set_att(0*dB)
        delay(1*ms)

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
        self.phaser.duc_stb()
        delay(1*ms)
        for _ in range(10000):
            t_stamp = self.get_frame_timestamp()
            at_mu(t_stamp + 100000)
            ch0.oscillator[0].set_frequency(7 * MHz)
            ch1.oscillator[0].set_frequency(7 * MHz)
            delay_mu(40*2)
            ch0.oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            ch1.oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            delay_mu(40*10)
            ch0.oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=0)
            ch1.oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=0)
            delay_mu(40*2)
            ch0.oscillator[0].set_amplitude_phase(amplitude=.8)
            ch1.oscillator[0].set_amplitude_phase(amplitude=.8)
            delay_mu(40*25)
            ch0.oscillator[0].set_amplitude_phase(amplitude=.8, phase=0.5)
            ch1.oscillator[0].set_amplitude_phase(amplitude=.8, phase=0.5)
            delay_mu(40*25)
            ch0.oscillator[0].set_amplitude_phase(amplitude=.0)
            ch1.oscillator[0].set_amplitude_phase(amplitude=.0)

            self.core.wait_until_mu(now_mu())
            print("make jitter")
            self.core.break_realtime()
            delay(10*ms)
