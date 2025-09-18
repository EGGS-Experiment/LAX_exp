from artiq.experiment import *
from artiq.coredevice import ad9910
from numpy import int32, int64, linspace

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence
from LAX_exp.system.objects.RAMWriter import RAMWriter
from LAX_exp.system.objects.PulseShaper import available_pulse_shapes

# Digital Ramp Generator - Ramp Destination
DRG_DEST_FTW =  0b00
DRG_DEST_POW =  0b01
DRG_DEST_ASF =  0b10

# todo: implement enable
# todo: also enable etc on qubitpulseshape etc.


class QubitRAP(LAXSubsequence):
    """
    Subsequence: Qubit Rapid Adiabatic Passage (RAP)

    Do rapid adiabatic passage via frequency-chirped + pulse-shaped
    coherent qubit pulse on 40Ca+ via the 729nm.
    """
    name = 'qubit_RAP'
    kernel_invariants = {
        # waveform configs
        "ram_profile", "ram_addr_start", "num_samples", "ampl_max_pct", "pulse_shape", "ram_addr_stop",

        # RAM-related configs
        "time_pulse_mu_to_ram_step", "time_pulse_mu_to_drg_step", "ampl_asf_pulseshape_list", "ram_writer",

        # register configurations
        "_CFR1_RAM_CONFIG", "_CFR2_DRG_CONFIG",
    }

    def build_subsequence(self, ram_profile: TInt32 = -1, ram_addr_start: TInt32 = 0x00,
                          num_samples: TInt32 = 200, ampl_max_pct: TFloat = 50.,
                          pulse_shape: TStr = "blackman", enable: TBool = True):
        """
        Defines the main interface for the subsequence.
        :param ram_profile: the AD9910 RAM profile to use for pulse shaping.
        :param ram_addr_start: the beginning RAM register address for pulse shaping.
            Must be in [0, 923].
        :param num_samples: the number of samples to use for the pulse shape.
            Must result in a final RAM address <= 1023.
        :param ampl_max_pct: the max amplitude (in percentage of full scale) of the pulse shape.
        :param pulse_shape: the pulse shape to use. Must be supported by available_pulse_shapes.
        :param enable: whether to actually enable the subsequence (via initialization etc.)
        """
        # set subsequence parameters
        self.ram_profile =      ram_profile
        self.ram_addr_start =   ram_addr_start
        self.num_samples =      num_samples
        self.ampl_max_pct =     ampl_max_pct
        self.pulse_shape =      pulse_shape
        self.enable =           enable

        # number of DRG updates per RAM amplitude update; must be power of 2
        self.drg_steps_per_ram_step = 1

        # get relevant devices
        self.setattr_device("core")
        self.setattr_device("qubit")
        self.ram_writer = RAMWriter(self, dds_device=self.qubit.beam,
                                    dds_profile=self.ram_profile, block_size=50)

    def prepare_subsequence(self):
        """
        Prepare values for speedy evaluation.
        """
        '''
        VALIDATE INPUTS
        '''
        self._prepare_argument_checks()


        '''SPECFIY RAM PARAMETERS FOR PULSE SHAPING'''
        # convert timings to multiples of SYNC_CLK (i.e. waveform update clock) period
        self.time_pulse_mu_to_ram_step = (self.qubit.beam.sysclk_per_mu / 4) / self.num_samples # SYNC_CLK period is 4x AD9910's SYSCLK
        self.time_pulse_mu_to_drg_step = self.time_pulse_mu_to_ram_step / self.drg_steps_per_ram_step

        # prepare ram_writer (b/c only LAXExperiment classes call their own children)
        self.ram_writer.prepare()
        self.ram_addr_stop = self.ram_addr_start + (self.num_samples - 1)

        # calculate pulse shape, then normalize and rescale to max amplitude
        # note: make sure max x_val is double the rolloff since PulseShaper does rising edge only
        x_vals = linspace(0., 200., self.num_samples)
        wav_y_vals = available_pulse_shapes[self.pulse_shape](x_vals, 100)
        wav_y_vals *= (self.ampl_max_pct / 100.) / max(wav_y_vals)

        # create array to store amplitude waveform in ASF (but formatted as a RAM word)
        self.ampl_asf_pulseshape_list = [int32(0)] * self.num_samples
        self.qubit.amplitude_to_ram(wav_y_vals, self.ampl_asf_pulseshape_list)
        self.ampl_asf_pulseshape_list = self.ampl_asf_pulseshape_list[::-1] # pre-reverse list b/c write_ram reverses it


        '''PREPARE CFR CONFIGURATION WORDS'''
        # CFR1: enable RAM mode and clear phase accumulator
        # note: has to be int64 b/c numpy won't take it as int32
        self._CFR1_RAM_CONFIG = int64(
            (1 << 31) | # ram_enable
            (ad9910.RAM_DEST_ASF << 29) |   # ram_destination
            (1 << 13) | # phase_autoclear
            (1 << 15) | # load_lrr
            (1 << 14) | # drg_autoclear
            2 # sdio_input_only + msb_first
        ) & 0xFFFFFFFF

        # CFR2: enable digital ramp generation
        # note: has to be int64 b/c numpy won't take it as int32
        # note: DRG nodwell low necessary for negative slopes b/c DRG
        #   accumulator always initialized to the lower limit
        self._CFR2_DRG_CONFIG = int64(
            (1 << 24) | # asf_profile_enable
            (1 << 16) | # effective_ftw
            (1 << 7) |  # matched_latency_enable
            (DRG_DEST_FTW << 20) |  # digital_ramp_destination
            (1 << 19) | # digital_ramp_enable
            (1 << 17) | # digital_ramp_nodwell_low
            (1 << 18)   # digital_ramp_nodwell_high
        ) & 0xFFFFFFFF

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # note: validate inputs here to get around bugs where args are passed from setattr_argument
        # sanitize sequence parameters
        # note: MUST USE PROFILE0 FOR BIDIRECTIONAL RAMP
        if self.ram_profile not in range(0, 7):
            raise ValueError("Invalid AD9910 profile for qubit_pulseshape: {:d}. Must be in [0, 7].".format(self.ram_profile))
        elif not (0 <= self.ram_addr_start <= 1023 - 100):
            raise ValueError("Invalid RAM start address for qubit_pulseshape: {:d}. Must be in [0, 923].".format(self.ram_addr_start))
        elif not (20 <= self.num_samples <= 1023 - self.ram_addr_start):
            raise ValueError("Invalid num_samples for qubit_pulseshape: {:d}. Must be in [20, 1000].".format(self.num_samples))
        elif not (0. <= self.ampl_max_pct <= 50.):
            raise ValueError("Invalid ampl_max_pct value ({:f}). Must be in range [0., 50.].".format(self.ampl_max_pct))


    """
    KERNEL FUNCTIONS
    """
    @kernel(flags={"fast-math"})
    def initialize_subsequence(self) -> TNone:
        """
        Prepare the subsequence immediately before run.
        """
        # disable RAM + DRG and set matched latencies
        self.qubit.set_cfr1(ram_enable=0)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.qubit.set_cfr2(matched_latency_enable=1)
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

        # write waveform to RAM via RAMWriter
        self.ram_writer.write(self.ampl_asf_pulseshape_list, self.ram_addr_start)
        delay_mu(25000)

        # set up RAM profile correctly after waveform uploaded
        self.qubit.set_profile_ram(
            start=self.ram_addr_start, end=self.ram_addr_stop,
            step=0xFFF,  # note: step size irrelevant since it's set in configure()
            profile=self.ram_profile, mode=ad9910.RAM_MODE_RAMPUP
        )
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(10000)

    @kernel(flags={"fast-math"})
    def cleanup_subsequence(self) -> TNone:
        """
        Clean up the subsequence immediately after run.
        """
        # stop & clear output
        self.qubit.off()
        self.qubit.set_ftw(0x00)
        self.qubit.set_asf(0x00)
        self.qubit.set_pow(0x00)
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

        # disable RAM mode
        self.qubit.set_cfr1(ram_enable=0)
        self.qubit.set_cfr2(matched_latency_enable=1)
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(256)   # add extra slack to avoid RTIO collisions

    @kernel(flags={"fast-math"})
    def configure(self, time_mu: TInt64, freq_center_ftw: TInt32, freq_dev_ftw: TInt32) -> TInt64:
        """
        Set the overall pulse time for the shaped pulse.
            This is achieved by adjusting the sample rate of the pulse shape updates.
        :param time_mu: pulse time in mu.
        :param freq_center_ftw: chirp center freq in 32b Frequency Tuning Word (ftw).
        :param freq_dev_ftw: chirp deviation (both directions) in 32b Frequency Tuning Word (ftw).
        :return: actual pulse time (in mu).
        """
        '''CALCULATE VALUES'''
        # calculate step size/timing for RAM
        time_ram_step = round(time_mu * self.time_pulse_mu_to_ram_step)
        if (time_ram_step > (1 << 16)) or (time_ram_step < 1):
            raise ValueError("Invalid RAM timestep in qubitRAP.configure. Change either pulse time or num_samples.")

        # calculate step size/timing for DRG
        time_drg_step = round(time_mu * self.time_pulse_mu_to_drg_step)
        if (time_drg_step > (1 << 16)) or (time_drg_step < 1):
            raise ValueError("Invalid DRG timestep in qubitRAP.configure. Change either pulse time or num_samples.")
        # note: freq_dev_ftw << 1 to make it double-sided
        freq_drg_ftw = int32(round((freq_dev_ftw << 1) / (self.num_samples * self.drg_steps_per_ram_step)))


        '''CONFIGURE HARDWARE - PULSE SHAPING'''
        # set RAM profile parameters for pulse shaping
        # note: using RAM rampup mode for simplicity
        self.qubit.set_profile_ram(
            start=self.ram_addr_start, end=self.ram_addr_stop,
            step=time_ram_step,
            profile=self.ram_profile, mode=ad9910.RAM_MODE_RAMPUP
        )
        self.qubit.cpld.io_update.pulse_mu(8)


        '''CONFIGURE HARDWARE - FREQUENCY RAMP/CHIRP'''
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

        # return relevant values
        return int64(time_ram_step / self.time_pulse_mu_to_ram_step)

    @kernel(flags={"fast-math"})
    def run_rap(self, time_pulse_mu: TInt64) -> TNone:
        """
        Fire pulse-shaped RAM pulse.
        :param time_pulse_mu: total pulse time.
            Doesn't have to be same as RAP pulse length (i.e. can cut off pulse early).
        """
        '''PRIME RAM + DRG'''
        # set target DDS profile
        self.qubit.cpld.set_profile(self.ram_profile)
        self.qubit.cpld.io_update.pulse_mu(8)

        # prepare CFRs
        self.qubit.write32(ad9910._AD9910_REG_CFR1, int32(self._CFR1_RAM_CONFIG)) # enable RAM & clear phase accumulator
        self.qubit.write32(ad9910._AD9910_REG_CFR2, int32(self._CFR2_DRG_CONFIG)) # enable digital ramp generation
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(256)   # necessary to prevent RTIO collisions

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

    @portable(flags={"fast-math"})
    def configure_values(self, time_mu: TInt64, freq_center_ftw: TInt32, freq_dev_ftw: TInt32
                      ) -> TTuple([TList(TInt32), TInt64]):
        """
        Precalculate register words for RAP (i.e. RAM profile, DRG limits/step/rate)
        :param time_mu: pulse time in mu.
        :param freq_center_ftw: chirp center freq in 32b Frequency Tuning Word (ftw).
        :param freq_dev_ftw: chirp deviation (both directions) in 32b Frequency Tuning Word (ftw).
        :return: list of register words, i.e. [ram_prof_hi, ram_prof_lo, drg_lim_hi, drg_lim_lo,
            drg_ramp_step_hi, drg_ramp_step_lo, drg_ramp_rate].
        :return: actual pulse time (in mu).
        """
        '''CALCULATE VALUES'''
        # calculate step size/timing for RAM
        time_ram_step = round(time_mu * self.time_pulse_mu_to_ram_step)
        if (time_ram_step > (1 << 16)) or (time_ram_step < 1):
            raise ValueError("Invalid RAM timestep in qubitRAP.configure. Change either pulse time or num_samples.")
        time_pulse_actual = int64(time_ram_step / self.time_pulse_mu_to_ram_step) # actual pulse time in hardware

        # calculate step size/timing for DRG
        time_drg_step = round(time_mu * self.time_pulse_mu_to_drg_step)
        if (time_drg_step > (1 << 16)) or (time_drg_step < 1):
            raise ValueError("Invalid DRG timestep in qubitRAP.configure. Change either pulse time or num_samples.")
        # note: freq_dev_ftw << 1 to make it double-sided
        freq_drg_ftw = int32(round((freq_dev_ftw << 1) / (self.num_samples * self.drg_steps_per_ram_step)))

        '''CONFIGURE REGISTER WORDS'''
        # RAM profile word - 64b: goes to _AD9910_REG_PROFILE0 + profile number
        # note: use RAM ramp up mode for simplicity
        ram_prof_hi = (time_ram_step << 8) | (self.ram_addr_stop >> 2)
        ram_prof_lo = ((self.ram_addr_stop << 30) | (self.ram_addr_start << 14) | (0 << 5) | (0 << 3) | ad9910.RAM_MODE_RAMPUP)

        # DRG ramp limits - 64b: goes to _AD9910_REG_RAMP_LIMIT
        drg_lim_hi = freq_center_ftw + freq_dev_ftw
        drg_lim_lo = freq_center_ftw - freq_dev_ftw

        # DRG ramp step size - 64b: goes to _AD9910_REG_RAMP_STEP
        drg_ramp_step_hi = -freq_drg_ftw
        drg_ramp_step_lo = freq_drg_ftw

        # DRG ramp rate - 32b: goes to _AD9910_REG_RAMP_RATE
        drg_ramp_rate = (time_drg_step << 16) | (time_drg_step << 0)

        # return relevant values
        return ([ram_prof_hi, ram_prof_lo, drg_lim_hi, drg_lim_lo, drg_ramp_step_hi, drg_ramp_step_lo, drg_ramp_rate],
                time_pulse_actual)

    @kernel(flags={"fast-math"})
    def run_from_config(self, register_words: TList(TInt32), time_pulse_mu: TInt64) -> TNone:
        """
        Fire pulse-shaped RAM pulse.
        :param register_words: list of 32b register words.
            Should be a length 7 list of [ram_prof_hi, ram_prof_lo, drg_lim_hi, drg_lim_lo,
            drg_ramp_step_hi, drg_ramp_step_lo, drg_ramp_rate].
        :param time_pulse_mu: total pulse time.
            Doesn't have to be same as RAP pulse length (i.e. can cut off pulse early).
        """
        # sanitize input (?)
        if len(register_words) != 7: raise ValueError("Invalid register_words. Should be a length 7 list.")

        '''CONFIGURE REGISTERS'''
        # set RAM profile parameters for amplitude pulse shaping
        self.qubit.write64(ad9910._AD9910_REG_PROFILE0 + self.ram_profile,
                           data_high=register_words[0],
                           data_low=register_words[1])
        self.qubit.cpld.io_update.pulse_mu(8)   # ensure profile is latched
        # set Digital Ramp Generator ramp limits
        self.qubit.write64(ad9910._AD9910_REG_RAMP_LIMIT,
                           data_high=register_words[2],  # max freq
                           data_low=register_words[3])   # min freq
        # set Digital Ramp Generator frequency step size
        self.qubit.write64(ad9910._AD9910_REG_RAMP_STEP,
                           data_high=register_words[4], # ramp down step size
                           data_low=register_words[5])  # ramp up step size
        # set Digital Ramp Generator update interval
        self.qubit.write32(ad9910._AD9910_REG_RAMP_RATE,
                           register_words[6])           # ramp up [16:8] & ramp down [7:0]
        self.qubit.cpld.io_update.pulse_mu(8)

        '''PRIME RAM + DRG'''
        # set target DDS profile
        self.qubit.cpld.set_profile(self.ram_profile)
        self.qubit.cpld.io_update.pulse_mu(8)

        # prepare CFRs
        self.qubit.write32(ad9910._AD9910_REG_CFR1, int32(self._CFR1_RAM_CONFIG)) # enable RAM & clear phase accumulator
        self.qubit.write32(ad9910._AD9910_REG_CFR2, int32(self._CFR2_DRG_CONFIG)) # enable digital ramp generation
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(256)   # necessary to prevent RTIO collisions

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

