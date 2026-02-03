import numpy as np
from artiq.experiment import *

from artiq.coredevice import ad9910, urukul
from LAX_exp.system.subsequences.dds_pulse_shape import DDSPulseShape


class DDSPulseshapeTest(EnvExperiment):
    """
    Test DDS pulse shaping.
    """
    kernel_invariants = {
        "profile_target", "ampl_pct", "freq_mhz",

        "dds_sig", "dds_ref", "ttl",

        "ftw", "asf", "att_db", "time_pulse_mu", "time_delay_mu",
    }

    def build(self):
        self.setattr_argument("repetitions", NumberValue(default=10000, precision=0, step=1, min=1, max=10000000))
        self.setattr_argument("time_pulse_us",  NumberValue(default=50, precision=3, step=500, min=1, max=10000000, unit='us', scale=1.))
        self.setattr_argument("time_delay_us",  NumberValue(default=5000, precision=3, step=500, min=1, max=10000000, unit='us', scale=1.))

        self.setattr_device("core")
        self.dds_sig =  self.get_device('urukul1_ch1')
        self.dds_ref =  self.get_device('urukul1_ch2')
        self.ttl =      self.get_device("ttl9")

        self.profile_target =   5
        self.ampl_pct =         50.
        self.freq_mhz =         120.339

    def prepare(self):
        self.ftw = self.dds.frequency_to_ftw(self.freq_mhz)
        self.asf = self.dds.amplitude_to_asf(self.ampl_pct / 100.)
        self.att_db = 7 * dB

        self.time_pulse_mu = self.core.seconds_to_mu(self.time_pulse_us * us)
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_us * us)
        self.pulseshape = DDSPulseShape(
            self,
            dds_target=self.dds_sig,
            ram_profile=self.profile_target,
            ram_addr_start=0,
            num_samples=500,
            ampl_max_pct=self.ampl_pct,
            pulse_shape="blackman",
        )


    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        self._run_initialize()

        # add slack
        self.core.break_realtime()
        delay_mu(1000000)

        # configure
        t_begin_mu = now_mu() & ~0x7
        self.dds_ref.set_mu(self.ftw, asf=self.asf, pow_=0x0,
                            profile=self.profile_target,
                            phase_mode=ad9910.PHASE_MODE_TRACKING,
                            ref_time_mu=t_begin_mu)
        self.dds_sig.set_ftw(self.ftw)
        self.pulseshape.configure_train(self.time_pulse_mu)

        # run pulse train
        for i in range(self.repetitions):
            self.ttl.on()
            self.pulseshape.run_train_single(phas_pow=0x0)
            self.ttl.off()
            delay_mu(self.time_delay_mu)


        self._run_cleanup()

    @kernel(flags={"fast-math"})
    def _run_initialize(self):
        self.core.reset()
        self.ttl.off()
        self.dds_sig.sw.off()
        self.dds_ref.sw.off()

        self.dds_sig.cpld.get_att_mu()
        self.core.break_realtime()
        self.dds_ref.cpld.get_att_mu()
        self.core.break_realtime()

        self.pulseshape.sequence_initialize()

        self.dds_sig.set_att(self.att_db)
        self.dds_ref.set_att(self.att_db)
        self.core.break_realtime()

        self.dds_sig.cpld.set_profile(self.profile_target)
        self.dds_ref.cpld.set_profile(self.profile_target)
        self.core.break_realtime()

        self.dds_sig.set_cfr2(matched_latency_enable=1)
        self.dds_ref.set_cfr2(matched_latency_enable=1)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def _run_cleanup(self):
        self.core.break_realtime()
        self.pulseshape.sequence_cleanup()

        self.ttl.off()
        self.dds_sig.sw.off()
        self.dds_ref.sw.off()

        self.dds_sig.set_att(31.5)
        self.dds_ref.set_att(31.5)
        self.core.break_realtime()


        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

