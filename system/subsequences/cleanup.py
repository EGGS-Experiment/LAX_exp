from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class Cleanup(LAXSubsequence):
    """
    Subsequence: Cleanup

    Reset device states to their default for client use after a sequence has been run.
    """
    name = 'cleanup'

    def build_subsequence(self):
        # get devices
        self.setattr_device('urukul0_cpld')
        self.setattr_device('urukul1_cpld')
        self.setattr_device('urukul2_cpld')
        self.setattr_device('phaser0')
        self.setattr_device('phaser1')
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')
        self.setattr_device('ttl10')
        self.setattr_device('ttl11')
        self.setattr_device('ttl12')
        self.setattr_device('ttl13')
        self.setattr_device('ttl14')

        # specific device cases
        self.setattr_device('urukul1_ch1')      # parametric DDS
        self.setattr_device('urukul1_ch2')      # tickle DDS

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        """
        Reset all ARTIQ hardware to ensure experiment is in a safe state.
        """
        # reset core device, RTIOs, and FIFOs
        self.core.reset()

        '''
        DDS HARDWARE
        '''
        # enable all RF switches
        self.ttl12.off()
        self.ttl13.off()
        self.ttl14.off()

        # reset qubit board
        self.urukul0_cpld.set_profile(0)
        self.urukul0_cpld.io_update.pulse_mu(8)
        self.urukul0_cpld.cfg_switches(0b0000)
        delay_mu(50000)
        # todo: ensure urukul0 TTL switches are also closed

        # reset motional board to rescue parameters
        self.urukul1_cpld.set_profile(0)
        self.urukul1_cpld.io_update.pulse_mu(8)
        self.urukul1_cpld.cfg_switches(0b0000)
        # set maximum attenuation for motional DDSs to prevent leakage
        self.urukul1_ch1.set_att_mu(0)
        self.urukul1_ch2.set_att_mu(0)
        # self.urukul1_cpld.set_all_att_mu(0)
        # set clean waveform for tickle DDS to prevent leakage
        # self.urukul1_ch3.set_mu(0x01, asf=0x01, profile=0)
        delay_mu(50000)
        # todo: ensure urukul1 TTL switches are also closed

        # reset main board to rescue parameters
        self.urukul2_cpld.set_profile(0)
        self.urukul2_cpld.io_update.pulse_mu(8)
        self.urukul2_cpld.cfg_switches(0b1110)
        delay_mu(50000)


        '''
        PARAMETRIC HARDWARE
        '''
        # stop parametric signal via external switch
        self.ttl11.off()


        '''
        PHASER HARDWARE
        '''
        # ensure INT HOLD on RF servo is deactivated
        self.ttl10.off()
        # ensure EGGS amp switches are closed
        self.ttl8.off()
        self.ttl9.off()
        delay_mu(50000)

        ### PHASER0 ###
        # reset phaser attenuators
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_att_mu(0x00)
        delay_mu(40)
        self.phaser0.channel[1].set_att_mu(0x00)
        delay_mu(100000)

        # reset phaser oscillators
        for i in range(5):
            # synchronize to frame
            delay_mu(125000)
            at_mu(self.phaser0.get_next_frame_mu())

            # clear oscillator frequencies
            with parallel:
                self.phaser0.channel[0].oscillator[i].set_frequency(0.)
                self.phaser0.channel[1].oscillator[i].set_frequency(0.)
                delay_mu(40)

            # clear oscillator amplitudes
            with parallel:
                self.phaser0.channel[0].oscillator[i].set_amplitude_phase(amplitude=0.)
                self.phaser0.channel[1].oscillator[i].set_amplitude_phase(amplitude=0.)
                delay_mu(40)


        ### PHASER1 ###
        # reset phaser attenuators
        delay_mu(125000)
        at_mu(self.phaser1.get_next_frame_mu())
        self.phaser1.channel[0].set_att_mu(0x00)
        delay_mu(40)
        self.phaser1.channel[1].set_att_mu(0x00)
        delay_mu(125000)

        # reset phaser oscillators
        for i in range(5):
            # synchronize to frame
            delay_mu(125000)
            at_mu(self.phaser1.get_next_frame_mu())

            # clear oscillator frequencies
            with parallel:
                self.phaser1.channel[0].oscillator[i].set_frequency(0.)
                self.phaser1.channel[1].oscillator[i].set_frequency(0.)
                delay_mu(40)

            # clear oscillator amplitudes
            with parallel:
                self.phaser1.channel[0].oscillator[i].set_amplitude_phase(amplitude=0.)
                self.phaser1.channel[1].oscillator[i].set_amplitude_phase(amplitude=0.)
                delay_mu(40)

        # add slack, then synchronize timeline to ensure events submit before resetting
        delay_mu(125000)
        self.core.wait_until_mu(now_mu())
        self.core.reset()
