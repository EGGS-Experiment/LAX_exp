import numpy as np
from artiq.experiment import *
from artiq.coredevice import ad9910

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence
from LAX_exp.system.objects.PulseShaper import available_pulse_shapes

# Digital Ramp Generator - Ramp Destination
DRG_DEST_FTW =  0b00
DRG_DEST_POW =  0b01
DRG_DEST_ASF =  0b10


class QubitRAP(LAXSubsequence):
    """
    Subsequence: Qubit Rapid Adiabatic Passage (RAP)

    Do rapid adiabatic passage via frequency-chirped + pulse-shaped
    coherent qubit pulse on 40Ca+ via the 729nm.
    """
    name = 'qubit_RAP'
    kernel_invariants = {
        "ram_profile", "ram_addr_start", "num_samples", "ampl_max_pct", "pulse_shape",
        "ram_addr_stop", "freq_dds_sync_clk_hz", "time_pulse_mu_to_ram_step", "time_pulse_mu_to_drg_step",
        "ampl_asf_pulseshape_list",
    }

    def build_subsequence(self, ram_profile: TInt32 = 0, ram_addr_start: TInt32 = 0x00,
                          num_samples: TInt32 = 1000, ampl_max_pct: TFloat = 50.,
                          pulse_shape: TStr = "blackman"):
        """
        Defines the main interface for the subsequence.
        Arguments:
            ram_profile: the AD9910 RAM profile to use for pulse shaping.
            ram_addr_start: the beginning RAM register address for pulse shaping.
                Must be in [0, 923].
            num_samples: the number of samples to use for the pulse shape.
                Must result in a final RAM address <= 1023.
            ampl_max_pct: the max amplitude (in percentage of full scale) of the pulse shape.
        """
        # set subsequence parameters
        self.ram_profile =      ram_profile
        self.ram_addr_start =   ram_addr_start
        # self.num_samples =      num_samples
        self.num_samples =      1000
        self.ampl_max_pct =     ampl_max_pct
        self.pulse_shape =      pulse_shape

        # number of DRG updates per RAM amplitude update; must be power of 2
        self.drg_steps_per_ram_step = 1

        # get relevant devices
        self.setattr_device("core")
        self.setattr_device("qubit")

    def prepare_subsequence(self):
        """
        Prepare values for speedy evaluation.
        """
        '''
        VALIDATE INPUTS
        '''
        # note: validate inputs here to get around bugs where args are passed from setattr_argument
        # sanitize sequence parameters
        # note: MUST USE PROFILE0 FOR BIDIRECTIONAL RAMP
        if self.ram_profile not in range(0, 7):
            raise ValueError("Invalid AD9910 profile for qubit_pulseshape: {:d}. Must be in [0, 7].".format(self.ram_profile))
        elif not (0 <= self.ram_addr_start <= 1023 - 100):
            raise ValueError("Invalid RAM start address for qubit_pulseshape: {:d}. Must be in [0, 923].".format(self.ram_addr_start))
        elif not (100 <= self.num_samples <= 1023 - self.ram_addr_start):
            raise ValueError("Invalid num_samples for qubit_pulseshape: {:d}. Must be in [100, 1000].".format(self.num_samples))
        elif not (0. <= self.ampl_max_pct <= 50.):
            raise ValueError("Invalid ampl_max_pct value ({:f}). Must be in range [0., 50.].".format(self.ampl_max_pct))

        '''SPECIFY TIMING'''
        # convert specified waveform sample rate to multiples of the SYNC_CLK (i.e. waveform update clock) period
        # todo: get sync_clk from ad9910 device instead
        self.freq_dds_sync_clk_hz = 1e9 / 4.    # SYNC_CLK HAS 4ns PERIOD
        self.time_pulse_mu_to_ram_step =    ((self.freq_dds_sync_clk_hz / self.core.seconds_to_mu(1)) /
                                             self.num_samples)
        self.time_pulse_mu_to_drg_step =    ((self.freq_dds_sync_clk_hz / self.core.seconds_to_mu(1)) /
                                             (self.num_samples * self.drg_steps_per_ram_step))

        '''SPECFIY RAM PARAMETERS FOR PULSE SHAPING'''
        # stop RAM address
        self.ram_addr_stop = self.ram_addr_start + (self.num_samples - 1)

        # calculate pulse shape, then normalize and rescale to max amplitude
        # note: make sure max x_val is double the rolloff since PulseShaper does rising edge only
        x_vals = np.linspace(0., 200., self.num_samples)
        wav_y_vals = available_pulse_shapes[self.pulse_shape](x_vals, 100)
        wav_y_vals *= (self.ampl_max_pct / 100.) / np.max(wav_y_vals)

        # create empty array to store values
        self.ampl_asf_pulseshape_list = [np.int32(0)] * self.num_samples
        # convert amplitude data to RAM in ampl. mod. mode (i.e. 64-bit word) and store in ampl_asf_pulseshape_list
        self.qubit.amplitude_to_ram(wav_y_vals, self.ampl_asf_pulseshape_list)
        # pre-reverse ampl_asf_pulseshape_list since write_ram makes a booboo and reverses the array
        self.ampl_asf_pulseshape_list = self.ampl_asf_pulseshape_list[::-1]


    """
    KERNEL FUNCTIONS
    """
    @kernel(flags={"fast-math"})
    def initialize_subsequence(self) -> TNone:
        """
        Prepare the subsequence immediately before run.
        """
        self.core.break_realtime()

        # disable RAM + DRG and set matched latencies
        self.qubit.set_cfr1(ram_enable=0)
        self.qubit.set_cfr2(matched_latency_enable=1)
        self.core.break_realtime()
        # note: somehow this delay is critical
        delay_mu(1000000)

        # prepare to write waveform to RAM profile
        self.qubit.set_profile_ram(
            start=self.ram_addr_start, end=self.ram_addr_stop,
            step=0xFFF, # note: step size irrelevant since it's set in configure()
            profile=self.ram_profile, mode=ad9910.RAM_MODE_RAMPUP
        )
        self.core.break_realtime()
        delay_mu(1000000)

        # set target RAM profile
        self.qubit.cpld.set_profile(self.ram_profile)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # write waveform to RAM profile
        self.core.break_realtime()
        delay_mu(20000000)   # 20 ms
        # note: this IO_UPDATE is necessary for slack reasons (cf the critical 1ms delay above)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.qubit.write_ram(self.ampl_asf_pulseshape_list)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def cleanup_subsequence(self) -> TNone:
        """
        Clean up the subsequence immediately after run.
        """
        self.core.break_realtime()

        # stop & clear output
        self.qubit.off()
        self.qubit.set_ftw(0x00)
        self.qubit.set_asf(0x00)
        self.qubit.set_pow(0x00)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # disable RAM mode
        self.qubit.set_cfr1(ram_enable=0)
        # tood: does set_cfr1 need its own io_update? or can I latch ALL with a single io_update?
        self.qubit.set_cfr2(matched_latency_enable=1)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def configure(self, time_mu: TInt64, freq_center_ftw: TInt32,
                  freq_dev_ftw: TInt32) -> TInt64:
        """
        Set the overall pulse time for the shaped pulse.
        This is achieved by adjusting the sample rate of the pulse shape updates.
        Arguments:
            time_mu: pulse time in mu.
            freq_center_ftw: chirp center freq in 32b Frequency Tuning Word (ftw).
            freq_dev_ftw: chirp deviation (both directions) in 32b Frequency Tuning Word (ftw).
        Returns:
            actual pulse time (in mu).
        """
        '''CALCULATE VALUES'''
        # calculate step size/timing for RAM
        time_ram_step = round(time_mu * self.time_pulse_mu_to_ram_step)
        if (time_ram_step > (1 << 16)) or (time_ram_step < 1):
            raise ValueError("Invalid RAM timestemp in qubitRAP.configure."
                             "Change either pulse time or number of samples.")

        # calculate step size/timing for DRG
        time_drg_step = round(time_mu * self.time_pulse_mu_to_drg_step)
        if (time_drg_step > (1 << 16)) or (time_drg_step < 1):
            raise ValueError("Invalid DRG timestep in qubitRAP.configure."
                             "Change either pulse time or number of samples.")
        # note: freq_dev_ftw << 1 to make it double-sided
        freq_drg_ftw = np.int32(round((freq_dev_ftw << 1) / (self.num_samples * self.drg_steps_per_ram_step)))

        '''CONFIGURE HARDWARE'''
        # set RAM profile parameters for pulse shaping
        # note: using RAM rampup mode for simplicity
        self.qubit.set_profile_ram(
            start=self.ram_addr_start, end=self.ram_addr_stop,
            step=time_ram_step,
            profile=self.ram_profile, mode=ad9910.RAM_MODE_RAMPUP
        )
        self.qubit.cpld.io_update.pulse_mu(8)

        # set Digital Ramp Generator limits
        self.qubit.write64(ad9910._AD9910_REG_RAMP_LIMIT,
                           data_high=freq_center_ftw+freq_dev_ftw,  # max freq
                           data_low=freq_center_ftw-freq_dev_ftw)   # min freq
        # set Digital Ramp Generator update interval
        self.qubit.write32(ad9910._AD9910_REG_RAMP_RATE,
                           (time_drg_step << 16) |   # ramp down
                           (time_drg_step << 0))     # ramp up
        # set Digital Ramp Generator frequency step size
        self.qubit.write64(ad9910._AD9910_REG_RAMP_STEP,
                           data_high=-freq_drg_ftw,  # ramp down
                           data_low=freq_drg_ftw)  # ramp up
        self.qubit.cpld.io_update.pulse_mu(8)

        # # tmp remove
        # print(time_pulse_mu, "us pulse time")
        # self.core.break_realtime()
        # delay_mu(1000000)
        # print(time_ram_step, " ram steps")
        # self.core.break_realtime()
        # delay_mu(1000000)
        # print(self.qubit.ftw_to_frequency(freq_center_ftw+freq_dev_ftw) / MHz, " MHz upper bound")
        # self.core.break_realtime()
        # delay_mu(1000000)
        # print(self.qubit.ftw_to_frequency(freq_center_ftw-freq_dev_ftw) / MHz, " MHz lower bound")
        # self.core.break_realtime()
        # delay_mu(1000000)
        # print(time_drg_step, " drg steps")
        # self.core.break_realtime()
        # delay_mu(1000000)
        # print(self.qubit.ftw_to_frequency(freq_drg_ftw) / kHz, " kHz per step (DRG)")
        # self.core.break_realtime()
        # delay_mu(1000000)
        # print((time_pulse_mu / (time_drg_step*4)) * (self.qubit.ftw_to_frequency(freq_drg_ftw) / kHz), " MHz")
        # self.core.break_realtime()
        # delay_mu(1000000)
        # # tmp remove

        # return relevant values
        return np.int64(time_ram_step / self.time_pulse_mu_to_ram_step)

    @kernel(flags={"fast-math"})
    def run_rap(self, time_pulse_mu: TInt64) -> TNone:
        """
        Fire pulse-shaped RAM pulse.
        Arguments:
            time_pulse_mu: the total pulse time. Doesn't have to be the
                same as the RAP pulse length (i.e. can cut off the pulse early).
        """
        '''PRIME RAM + DRG'''
        # set target DDS profile
        self.qubit.cpld.set_profile(self.ram_profile)
        self.qubit.cpld.io_update.pulse_mu(8)

        # enable RAM mode and clear DDS phase accumulator
        self.qubit.write32(ad9910._AD9910_REG_CFR1,
                           (1 << 31) |  # ram_enable
                           (ad9910.RAM_DEST_ASF << 29) |   # ram_destination
                           (1 << 13) |  # phase_autoclear
                           (1 << 15) |  # load_lrr
                           (1 << 14)    # drg_autoclear
                        )
        # enable digital ramp generation
        # note: DRG nodwell low is necessary to allow negative slopes
        # since DRG accumulator is always initialized to the lower limit
        self.qubit.write32(ad9910._AD9910_REG_CFR2,
                         (1 << 24) |  # asf_profile_enable
                         (1 << 16) |  # effective_ftw
                         (1 << 7) |  # matched_latency_enable
                         (DRG_DEST_FTW << 20) |  # digital_ramp_destination
                         (1 << 19) |  # digital_ramp_enable
                         (1 << 17) |  # digital_ramp_nodwell_low
                         (1 << 18)  # digital_ramp_nodwell_high
                         )
        self.qubit.cpld.io_update.pulse_mu(8)

        '''FIRE PULSE'''
        time_start_mu = now_mu() & ~7

        # start ramp-up when coarse aligned to SYNC_CLK for determinacy
        at_mu(time_start_mu)
        self.qubit.cpld.io_update.pulse_mu(8)

        # open and close switch to synchronize with RAM pulse
        at_mu(time_start_mu + 416 + 63 - 140 - 244)
        self.qubit.on()
        delay_mu(time_pulse_mu)
        self.qubit.off()

        '''CLEANUP'''
        # disable RAM and DRG
        self.qubit.set_cfr1(ram_enable=0)
        self.qubit.set_cfr2(matched_latency_enable=1)
        self.qubit.cpld.io_update.pulse_mu(8)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        pass
