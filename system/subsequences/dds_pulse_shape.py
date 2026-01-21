from artiq.experiment import *
from artiq.coredevice import ad9910
from numpy import int32, int64, linspace

from LAX_exp.system.objects.RAMWriter import RAMWriter
from LAX_exp.system.objects.PulseShaper import available_pulse_shapes


class DDSPulseShape(HasEnvironment):
    """
    Sequence: DDS Pulse Shape

    Apply a coherent, shaped pulse via RAM modulation on the AD9910.
    """
    name = 'DDSPulseShape'
    kernel_invariants = {
        # objects
        "dds_target", "ram_writer",

        # waveform configs
        "ram_profile", "ram_addr_start", "num_samples", "ampl_max_pct", "pulse_shape", "_CFR1_RAM_CONFIG",

        # RAM-related configs
        "ram_addr_stop", "time_pulse_mu_to_ram_step", "ampl_asf_pulseshape_list",
    }

    def build(self, dds_target=None,
              ram_profile: TInt32=0, ram_addr_start: TInt32=0x00,
              num_samples: TInt32=100, ampl_max_pct: TFloat=50.,
              pulse_shape: TStr="sine_squared"):
        """
        Defines the main interface for the sequence.
        :param dds_target: the DDS object to pulse shape with.
        :param ram_profile: the AD9910 RAM profile to use for pulse shaping.
        :param ram_addr_start: the beginning RAM register address for pulse shaping. Must be in [0, 923].
        :param num_samples: the number of samples to use for the pulse shape.
            Must result in a final RAM address <= 1023.
        :param ampl_max_pct: the max amplitude (in percentage of full scale) of the pulse shape.
        :param pulse_shape: the pulse shape to use.

        """
        # set sequence parameters
        self.ram_profile =      ram_profile
        self.ram_addr_start =   ram_addr_start
        self.num_samples =      num_samples
        self.ampl_max_pct =     ampl_max_pct
        self.pulse_shape =      pulse_shape

        # get relevant devices
        self.dds_target = dds_target
        self.setattr_device("core")
        self.ram_writer = RAMWriter(self, dds_device=self.dds_target,
                                    dds_profile=self.ram_profile, block_size=50)

        ### finish build sequence ###
        self._validate_arguments()
        self.sequence_prepare()


    """
    BUILD SEQUENCE
    """
    def _validate_arguments(self) -> TNone:
        """
        Check sequence arguments for validity.
        """
        # ensure dds_target is an AD9910 device
        if not isinstance(self.dds_target, ad9910.AD9910):
            raise ValueError("Invalid type for dds_target ({}). Must be of type AD9910.".format(
                type(self.dds_target)))

        # check parameters etc
        if self.ram_profile not in range(0, 7):
            raise ValueError("Invalid AD9910 profile for DDSPulseShape: {:d}. Must be in [0, 7].".format(self.ram_profile))
        elif not (0 <= self.ram_addr_start <= (1023 - 100)):
            raise ValueError("Invalid RAM start address for DDSPulseShape: {:d}. Must be in [0, 923].".format(self.ram_addr_start))
        elif not (20 <= self.num_samples <= (1023 - self.ram_addr_start)):
            raise ValueError("Invalid num_samples for DDSPulseShape: {:d}. Must be in [20, 1000].".format(self.num_samples))
        elif not (0. < self.ampl_max_pct < 100.):
            raise ValueError("Invalid ampl_max_pct value ({:f}). Must be in range (0., 100.).".format(self.ampl_max_pct))

    def sequence_prepare(self):
        """
        Prepare & pre-compute relevant values for speedy compilation & kernel evaluation.
        """
        # prepare stuff to write to RAM as well as our personal RAMWriter object
        self.ram_writer.prepare()
        self.ram_addr_stop = self.ram_addr_start + (self.num_samples - 1)
        self.time_pulse_mu = int64(0)    # holder variable for delay time (for later use)

        # convert timings to multiples of SYNC_CLK (i.e. waveform update clock) period
        self.time_pulse_mu_to_ram_step = ((self.dds_target.sysclk_per_mu / 4) /
                                          self.num_samples) # SYNC_CLK period is 4x AD9910's SYSCLK


        '''
        CALCULATE WAVEFORM
        '''
        # calculate pulse shape, then normalize and rescale to max amplitude
        # note: ensure max x_val is double the rolloff since PulseShaper does rising edge only
        x_vals = linspace(0., 200., self.num_samples)
        wav_y_vals = available_pulse_shapes[self.pulse_shape](x_vals, 100)
        wav_y_vals *= (self.ampl_max_pct / 100.) / max(wav_y_vals)

        # create array to store amplitude waveform in ASF (but formatted as a RAM word)
        self.ampl_asf_pulseshape_list = [int32(0)] * self.num_samples
        self.dds_target.amplitude_to_ram(wav_y_vals, self.ampl_asf_pulseshape_list)
        self.ampl_asf_pulseshape_list = self.ampl_asf_pulseshape_list[::-1] # pre-reverse list b/c write_ram reverses it


        '''
        PREPARE CFR CONFIGURATION WORDS
        '''
        # CFR1: enable RAM mode and clear phase accumulator
        # note: has to be int64 b/c numpy won't take it as int32
        self._CFR1_RAM_CONFIG = int64(
            (1 << 31) | # ram_enable
            (ad9910.RAM_DEST_ASF << 29) |   # ram_destination
            # (1 << 16) |  # select_sine_output (note: removed for consistency)
            2 # sdio_input_only + msb_first
        ) & 0xFFFFFFFF # ensure 32b only


    """
    INITIALIZE & CLEAN UP
    """
    @kernel(flags={"fast-math"})
    def sequence_initialize(self) -> TNone:
        """
        Prepare experiment hardware for this sequence.
        Should be called when the parent experiment enters the kernel.
        """
        self.core.break_realtime()

        # disable RAM mode & set matched latencies
        self.dds_target.set_cfr1(ram_enable=0)
        self.dds_target.cpld.io_update.pulse_mu(8)
        self.dds_target.set_cfr2(matched_latency_enable=1)
        self.dds_target.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

        # write waveform to RAM via RAMWriter
        self.ram_writer.write(self.ampl_asf_pulseshape_list, self.ram_addr_start)
        delay_mu(25000)

        # set up RAM profile correctly after waveform uploaded
        self.dds_target.set_profile_ram(
            start=self.ram_addr_start, end=self.ram_addr_stop,
            step=0xFFF,
            profile=self.ram_profile, mode=ad9910.RAM_MODE_RAMPUP
        )
        self.dds_target.cpld.io_update.pulse_mu(8)
        delay_mu(10000)

        # clear ASF and POW registers
        self.dds_target.set_asf(0x00)
        self.dds_target.set_pow(0x00)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def sequence_cleanup(self) -> TNone:
        """
        Leave experiment hardware (specific to this sequence) in a safe state.
        Should be called when the parent experiment completes the main experiment sequence and finishes.
        """
        self.core.break_realtime()

        # stop output & clear registers
        self.dds_target.sw.off()
        self.dds_target.set_asf(0x00)
        self.dds_target.set_pow(0x00)
        self.dds_target.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

        # disable RAM mode & return CFR1 to normal
        self.dds_target.set_cfr1(ram_enable=0)
        self.dds_target.cpld.io_update.pulse_mu(8)
        delay_mu(256)   # add extra slack to avoid RTIO collisions
        self.core.break_realtime()


    """
    CONFIGURABLE METHODS
    """
    @portable(flags={"fast-math"})
    def calculate_time(self, time_mu: TInt64) -> TTuple([TInt32, TInt64]):
        """
        Precalculate timing words for a RAM by adjusting the sample rate of the RAM waveform updates.
        :param time_mu: target pulse time in mu.
        :return: tuple of (time_ram_step (32b int truncated to 16b), time_pulse_actual_mu (64b int))
        """
        # calculate step size/timing for RAM
        time_ram_step = round(time_mu * self.time_pulse_mu_to_ram_step)
        if (time_ram_step > (1 << 16)) or (time_ram_step < 1):
            raise ValueError("Invalid RAM timestep in DDSPulseShape.calculate_time.\n"
                             "Change either pulse time or num_samples.")
        time_pulse_actual = int64(time_ram_step / self.time_pulse_mu_to_ram_step) # actual pulse time in hardware

        return time_ram_step, time_pulse_actual

    @kernel(flags={"fast-math"})
    def run(self, time_ram_step: TInt32, time_pulse_mu: TInt64, phas_pow: TInt32=0) -> TNone:
        """
        # todo: note phase discontinuity
        Fire shaped RAM pulse using values passed in as arguments.
        This function does not use any precomputed instance variables and makes no assumptions
            about the hardware state (i.e. it sets everything up from scratch every time it runs).
        :param time_ram_step: the timesteps for each RAM sample
        :param time_pulse_mu: total pulse time.
            Doesn't have to be same as target pulse length (i.e. can stop pulse early).
        :param phas_pow: DDS phase (in pow).
        """
        ##### CONFIGURE HARDWARE #####
        # set target phase
        self.dds_target.set_pow(phas_pow)

        # set RAM profile parameters for amplitude pulse shaping
        self.dds_target.set_profile_ram(start=self.ram_addr_start, end=self.ram_addr_stop,
                                        profile=self.ram_profile, mode=ad9910.RAM_MODE_RAMPUP,
                                        step=time_ram_step)
        self.dds_target.cpld.io_update.pulse_mu(8)   # ensure profile is latched

        # set target DDS profile
        self.dds_target.cpld.set_profile(self.ram_profile)
        self.dds_target.cpld.io_update.pulse_mu(8)

        # prepare CFRs
        self.dds_target.write32(ad9910._AD9910_REG_CFR1, int32(self._CFR1_RAM_CONFIG)) # enable RAM
        self.dds_target.cpld.io_update.pulse_mu(8)
        delay_mu(256)   # necessary to prevent RTIO collisions


        ##### FIRE PULSE #####
        time_start_mu = now_mu() & ~7 # coarse align to SYNC_CLK for determinacy
        at_mu(time_start_mu)
        self.dds_target.cpld.io_update.pulse_mu(8)

        # open and close switch to synchronize with RAM pulse
        at_mu(time_start_mu + 416 + 63 - 140 - 244)
        self.dds_target.sw.on()
        delay_mu(time_pulse_mu)
        self.dds_target.sw.off()


        ##### CLEANUP #####
        # disable RAM and DRG
        self.dds_target.set_cfr1(ram_enable=0)
        self.dds_target.cpld.io_update.pulse_mu(8)


    """
    TRAIN METHODS
    """
    # todo: note rigor of phase continuity, assumptions made, focus on speed
    @kernel(flags={"fast-math"})
    def configure_train(self, time_mu: TInt64) -> TInt64:
        """
        Complete configuration/setup of the AD9910 for RAM mode.
        This function is designed for fast, coherent, pulse train operation together with run_train_single.
        Once this function is called, the AD9910 hardware will be left ready to run a bunch of pulse trains
            (i.e. once you call this function, you should call run_train_single very very soon).
            Hardware is left in RAM mode (though with no output) following this function to reduce setup overhead,
            and IMPORTANTLY, to preserve phase coherence.
        Relevant values are precalculated and stored as an instance variable for fast later operation.
        :param time_mu: desired pulse time (in mu) for run_train_single.
        :return: actual pulse time (in mu) that run_train_single will take.
        """
        # calculate step size/timing for RAM
        time_ram_step = round(time_mu * self.time_pulse_mu_to_ram_step)
        if (time_ram_step > (1 << 16)) or (time_ram_step < 1):
            raise ValueError("Invalid RAM timestep in DDSPulseShape.configure_train().\n"
                             "Change either pulse time or num_samples.")
        # reconvert to get correct time_pulse_mu correctly for later delay
        self.time_pulse_mu = int64(time_ram_step / self.time_pulse_mu_to_ram_step)

        # set RAM profile parameters for pulse shaping
        # note: using RAM rampup mode for simplicity
        self.dds_target.set_profile_ram(
            start=self.ram_addr_start, end=self.ram_addr_stop,
            step=time_ram_step,
            profile=self.ram_profile, mode=ad9910.RAM_MODE_RAMPUP
        )
        self.dds_target.cpld.io_update.pulse_mu(8)   # ensure profile is latched

        # switch to target profile
        self.dds_target.cpld.set_profile(self.ram_profile)
        self.dds_target.cpld.io_update.pulse_mu(8)   # ensure profile is latched

        # enable RAM mode and leave primed
        self.dds_target.write32(ad9910._AD9910_REG_CFR1, int32(self._CFR1_RAM_CONFIG))
        self.dds_target.cpld.io_update.pulse_mu(8)  # ensure profile is latched

        return self.time_pulse_mu

    @kernel(flags={"fast-math"})
    def run_train_single(self, phas_pow: TInt32=-1) -> TNone:
        """
        Fire shaped pulse via RAM modulation functionality.
        This function uses timing values computed by configure_fixed() and assumes that
            the target RAM Profile control register is correct (which it should be if the
            register hasn't been touched since configure_fixed() was last run).
        :param phas_pow: DDS phase (in pow).
        """
        # note: in the name of speed, we have the option to NOT write the pow register,
        #   which means that run times will be different depending on whether we pass this
        #   function an argument or not.
        if phas_pow != -1: self.dds_target.set_pow(phas_pow) # set target phase

        time_start_mu = now_mu() & ~7 # coarse align to SYNC_CLK for determinacy
        at_mu(time_start_mu)
        self.dds_target.cpld.io_update.pulse_mu(8) # fire pulse!

        # open and close switch to synchronize with RAM pulse
        at_mu(time_start_mu + 416 + 63 - 140 - 244)
        self.dds_target.sw.on()
        delay_mu(self.time_pulse_mu)
        self.dds_target.sw.off()

        # note: we don't clean up (i.e. unset CFR1 or change profiles) b/c we want the
        #   phase accumulator to keep counting - this means that subsequent pulses will
        #   be phase-coherent (i.e. have a deterministic and well-defined phase relationship)
        #   to this pulse.
