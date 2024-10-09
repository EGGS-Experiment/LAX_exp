import numpy as np
from artiq.experiment import *


class PositionSwitch(EnvExperiment):
    """
    Position Switch
    """
    kernel_invariants = {
        "dds", "dds_att", "dds_ampl", "dds_freq", "dds_profile", "dds_time", "dds_off_att",
        "pmt", "pmt_retrievals", "pmt_time",
        "time_recooling_mu", "time_flipper_trigger_mu", "time_flipper_wait_mu"
    }


    def build(self):
        # get core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        self.setattr_device("ttl0_counter")
        self.setattr_device("ttl15")
        self.setattr_device("urukul1_ch3")


    def prepare(self):
        self.dds =          self.get_device("urukul1_ch3")
        self.dds_att =      self.dds.cpld.att_to_mu(0. * dB)
        self.dds_ampl =     self.dds.amplitude_to_asf(0.98)
        self.dds_freq =     self.dds.frequency_to_ftw(1373.5 * kHz)
        self.dds_profile =  4

        self.dds_time =     self.core.seconds_to_mu(10 * ms)
        self.dds_off_att =  self.dds.cpld.att_to_mu(31.5 * dB)

        # self._th0 = self.core.precompile(self.shuffle)
        # self._th1 = self.core.precompile(self.cleanup)

        self.pmt =              self.ttl0_counter
        self.pmt_time =         self.core.seconds_to_mu(50000 * us)
        self.pmt_retrievals =   10
        self.pmt_shuffle_threshold = 3100
        self._pmt_count_store = 0

        self.time_recooling_mu =        self.core.seconds_to_mu(1000 * ms)
        self.time_flipper_trigger_mu =  self.core.seconds_to_mu(50 * ms)
        self.time_flipper_wait_mu =     self.core.seconds_to_mu(2 * s)

    @kernel(flags={"fast-math"})
    def run(self):
        self.run_prepare()
        # todo: only run for x reps

        # check counts
        res0 = self.pmt_read(self.pmt_retrievals)
        self.core.break_realtime()
        print(res0)
        self.core.break_realtime()
        self.core.break_realtime()

        # main loop
        _loop_counter = 0
        while (res0 < self.pmt_shuffle_threshold) and (_loop_counter < 10):

            # shuffle
            self.shuffle()
            self.core.break_realtime()

            # add recooling time
            delay_mu(self.time_recooling_mu)

            # check counts
            res0 = self.pmt_read(self.pmt_retrievals)
            self.core.break_realtime()
            print(res0)
            self.core.break_realtime()
            self.core.break_realtime()

            # increment loop counter
            _loop_counter += 1

        # clean up
        self.core.break_realtime()
        if _loop_counter < 10:
            print('\n\tDESIRED POSITION!')
        else:
            print('\n\tSECULAR FREQUENCY CHANGE.')
        self.core.break_realtime()
        self.cleanup()

    @kernel(flags={"fast-math"})
    def run_prepare(self):
        self.core.break_realtime()
        self.core.wait_until_mu(now_mu())
        self.core.reset()

        # check whether we need to flip to PMT
        base_counts = self.pmt_read(10)
        self.core.break_realtime()
        if base_counts < 50:
            self.pmt_flip()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def pmt_read(self, num_retrievals: TInt32) -> TFloat:
        self.core.break_realtime()

        # reset PMT count holder
        self._pmt_count_store = 0

        # get PMT counts
        for i in range(num_retrievals):
            self.pmt.gate_rising_mu(self.pmt_time)
            self._pmt_count_store += self.pmt.fetch_count()
            delay_mu(25000)
        self.core.break_realtime()

        return self._pmt_count_store / num_retrievals


    @kernel(flags={"fast-math"})
    def pmt_flip(self):
        # flip to PMT
        self.ttl15.on()
        delay_mu(self.time_flipper_trigger_mu)
        self.ttl15.off()
        delay_mu(self.time_flipper_wait_mu)


    @kernel(flags={"fast-math"})
    def shuffle(self):
        self.core.break_realtime()

        # set waveform
        self.dds.set_att_mu(self.dds_att)
        self.dds.set_mu(self.dds_freq, asf=self.dds_ampl, profile=self.dds_profile)
        self.core.break_realtime()

        # shuffle
        self.dds.cpld.set_profile(self.dds_profile)
        self.dds.sw.on()
        delay_mu(self.dds_time)
        self.dds.sw.off()

        # clean up to prevent DDS overheating
        self.dds.set_att_mu(self.dds_off_att)
        self.dds.set_mu(self.dds_freq, asf=0x01, profile=self.dds_profile)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def cleanup(self):
        self.core.break_realtime()

        # set urukul
        self.dds.sw.off()
        self.dds.set_att_mu(self.dds_off_att)
        self.dds.set_mu(self.dds_freq, asf=0x01, profile=self.dds_profile)

        # flip to camera
        self.pmt_flip()

        # wait until end
        self.core.break_realtime()
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()
