from cmath import isnan

from artiq.experiment import *
from artiq.coredevice import ad9910, ttl
from numpy import int32, int64, linspace

from LAX_exp.system.objects.RAMWriter import RAMWriter
from LAX_exp.system.objects.PulseShaper import available_pulse_shapes
from LAX_exp.base import LAXDevice


class DDSPulseShaper(HasEnvironment):
    """
    Sequence: DDS Pulse Shaper

    Apply a coherent, shaped pulse via RAM modulation on the AD9910.
    """
    name = 'DDSPulseShaper'
    kernel_invariants = {
        # objects
        "dds_targets", "ram_writers",

        # waveform configs
        "ram_profile", "ram_addr_start", "num_samples", "ampl_max_pcts", "pulse_shapes",

        # RAM-related configs
        "ram_addr_stop", "time_pulse_mu_to_ram_step", "ampl_asf_pulseshapes_list",
    }

    def build(self, dds_target: ad9910.AD9910  = None,
              ram_profile: TInt32=0, ram_addr_start: TInt32=0x00,
              num_samples: TInt32=100, ampl_max_pct:  TFloat =50.,
              pulse_shape:  TStr= "sine_squared",
              pulse_configuration = 'all',
              phase_autoclear:  TInt32 = 0,
              external_switch: ttl.TTLOut = None):
        """
        Defines interface for pulse shaping multiple ddses on a single urukul board (share a CPLD)
        :param dds_targets: DDS objects to pulse shape with.
        :param ram_profile: the AD9910 RAM profile to use for pulse shaping (shared among DDSes)
        :param ram_addr_start: the beginning RAM register address for pulse shaping. Must be in [0, 923].
        :param num_samples: the number of samples to use for the pulse shape.
            Must result in a final RAM address <= 1023.
        :param ampl_max_pct: the max amplitude (in percentage of full scale) of the pulse shape.
        :param pulse_shape: the pulse shape to use.
        :param phase_autoclear: bit to determine whether to clear phase accumulator on io_update

        """


        # set sequence parameters
        self.ram_profile =      ram_profile
        self.ram_addr_start =   ram_addr_start
        self.num_samples =      num_samples
        self.ampl_max_pcts =     ampl_max_pct
        self.pulse_shapes =      pulse_shape
        self.external_switches = external_switch
        self.pulse_configuration = pulse_configuration

        # get relevant devices
        self.dds_targets = dds_target
        self.setattr_device("core")

        self._validate_arguments()


        self.ram_writers = [0] * len(self.dds_targets)
        for dds_targets_idx in range(len(self.dds_targets)):
            self.ram_writers[dds_targets_idx] = RAMWriter(self, dds_device=self.dds_targets[dds_targets_idx],
                                    dds_profile=self.ram_profile, block_size=50)

        self._cfr1_ram_configs = [0] * len(self.dds_targets)
        self.ampl_asf_pulseshapes_list = [[int32(0)] * self.num_samples]
        self.ampl_asf_pulseshapes_list_temp = [[int32(0)] * self.num_samples]


        ### finish build sequence ###
        self.sequence_prepare(phase_autoclear, 0)


    # add auxilary dds target
    def add_dds_target(self,
                       dds_target: ad9910.AD9910 = None,
                       ampl_max_pct: TFloat = 50.,
                       pulse_shape: TStr = "sine_squared",
                       phase_autoclear: TInt32 = 0,
                       external_switch: ttl.TTLOut = None):

        self.dds_targets +=  [dds_target]
        self.ampl_max_pcts += [ampl_max_pct]
        self.pulse_shapes += [pulse_shape]
        self.external_switches += [external_switch]
        self.ampl_asf_pulseshapes_list += [[int32(0)] * self.num_samples]
        self._validate_arguments()



        self.ram_writers += [RAMWriter(self, dds_device= dds_target,
                                                          dds_profile=self.ram_profile, block_size=50)]

        ### finish build sequence ###
        self.sequence_prepare(phase_autoclear,  len(self.dds_targets)-1)



    """
    BUILD SEQUENCE
    """
    def _validate_arguments(self) -> TNone:
        """
        Check sequence arguments for validity.
        """

        # if not in list format wrap in list
        if not isinstance(self.dds_targets, list) and (isinstance(self.dds_targets, ad9910.AD9910)
                or isinstance(self.dds_targets, LAXDevice)):
            self.dds_targets = [self.dds_targets]

        if not isinstance(self.ampl_max_pcts, list) and (isinstance(self.ampl_max_pcts, float)
                or isinstance(self.ampl_max_pcts, int32)):
            self.ampl_max_pcts = [self.ampl_max_pcts]

        if not isinstance(self.pulse_shapes, list) and isinstance(self.pulse_shapes, str):
                self.pulse_shapes = [self.pulse_shapes]

        if not isinstance(self.external_switches, list):
            if self.external_switches is None:
                self.external_switches = [0] * len(self.dds_targets)
            elif isinstance(self.external_switches, ttl.TTLOut):
                self.external_switches = [self.external_switches]

        # if only one pulse shape given and there are multiple dds targets assume we want to use the same pulse shape for all dds targets
        if len(self.pulse_shapes) == 1 and len(self.dds_targets) > 1:
            self.pulse_shapes = [self.pulse_shapes] * len(self.dds_targets)

        # if only one ampl max pct given and there are multiple dds targets assume we want to use the same max amplitude
        if len(self.ampl_max_pcts) == 1 and len(self.dds_targets) > 1:
            self.ampl_max_pcts = [self.ampl_max_pcts] * len(self.ampl_max_pcts)


        if not (len(self.dds_targets) == len(self.ampl_max_pcts) == len(self.pulse_shapes) == len(self.external_switches)):
            raise ValueError("All parameters must have the same shape or be single-valued")

        # ensure dds_target is an AD9910 device
        for dds_target in self.dds_targets:
            if not (isinstance(dds_target, ad9910.AD9910) or isinstance(dds_target, LAXDevice)):
                raise ValueError("Invalid type for dds_target ({}). Must be of type AD9910.".format(
                    type(dds_target)))

        # check parameters etc
        if self.ram_profile not in range(0, 7):
            raise ValueError("Invalid AD9910 profile for DDSPulseShape: {:d}. Must be in [0, 7].".format(self.ram_profile))
        elif not (0 <= self.ram_addr_start <= (1023 - 100)):
            raise ValueError("Invalid RAM start address for DDSPulseShape: {:d}. Must be in [0, 923].".format(self.ram_addr_start))
        elif not (20 <= self.num_samples <= (1023 - self.ram_addr_start)):
            raise ValueError("Invalid num_samples for DDSPulseShape: {:d}. Must be in [20, 1000].".format(self.num_samples))

        for ampl_max_pct in self.ampl_max_pcts:
            if not (0. < ampl_max_pct < 100.):
                raise ValueError("Invalid ampl_max_pct value ({:f}). Must be in range (0., 100.).".format(ampl_max_pct))

        for switch_idx in range(len(self.external_switches)):
            if self.external_switches[switch_idx] == 0 or self.external_switches[switch_idx] is None:
                self.external_switches[switch_idx] = self.DummyTTL()
            elif isinstance(self.external_switches[switch_idx], self.DummyTTL):
                pass
            elif not isinstance(self.external_switches[switch_idx], ttl.TTLOut):
                raise TypeError("the specified the external_switch must be an instance of TTLOut.")
            
        if self.pulse_configuration not in ['all', 'rising', 'falling']:
            raise ValueError('Need to configure the pulse properly')


    def sequence_prepare(self, phase_autoclear: TInt32, dds_targets_idx: TInt32) -> TNone:
        """
        Prepare & pre-compute relevant values for speedy compilation & kernel evaluation.

        :param phase_autoclear: whether to clear phase accumulator on io_update
        """
        # prepare stuff to write to RAM as well as our personal RAMWriter object
        self.time_pulse_mu_list = [int64(0)] * len(self.dds_targets)  # holder variable for delay time (for later use)
        dds_target = self.dds_targets[dds_targets_idx]
        self.ram_writers[dds_targets_idx].prepare()
        self.ram_addr_stop = self.ram_addr_start + (self.num_samples - 1)

        # convert timings to multiples of SYNC_CLK (i.e. waveform update clock) period
        self.time_pulse_mu_to_ram_step = ((dds_target.sysclk_per_mu / 4) /
                                          self.num_samples) # SYNC_CLK period is 4x AD9910's SYSCLK

        self.ram_firing_delay = int64(95)
        self.bit_shift_ram_enable = 31
        self.bit_shift_ram_dest = 29
        self.bit_shift_phase_autoclear = 13


        '''
        CALCULATE WAVEFORM
        '''
        # calculate pulse shape, then normalize and rescale to max amplitude
        # note: ensure max x_val is double the rolloff since PulseShaper does rising edge only
        x_vals = linspace(0., 200., self.num_samples)
        
        if self.pulse_configuration == 'all':
            samples_roll = int32(self.num_samples // 2)
        elif self.pulse_configuration == 'rising' or 'falling':
            samples_roll = int32(self.num_samples)
            
        wav_y_vals = available_pulse_shapes[self.pulse_shapes[dds_targets_idx]](x_vals, samples_roll)
        wav_y_vals *= (self.ampl_max_pcts[dds_targets_idx] / samples_roll) / max(wav_y_vals)

        # create array to store amplitude waveform in ASF (but formatted as a RAM word)
        self.ampl_asf_pulseshapes_list_temp = [int32(0)] * self.num_samples
        dds_target.amplitude_to_ram(wav_y_vals, self.ampl_asf_pulseshapes_list_temp)
        if self.pulse_configuration == 'rising' or self.pulse_configuration == 'all':
            self.ampl_asf_pulseshapes_list[dds_targets_idx] = self.ampl_asf_pulseshapes_list_temp # pre-reverse list b/c write_ram reverses it
        else:
            self.ampl_asf_pulseshapes_list[dds_targets_idx] = self.ampl_asf_pulseshapes_list_temp[::-1]  # pre-reverse list b/c write_ram reverses it

        # CFR1: enable RAM mode and clear phase accumulator
        # note: has to be int64 b/c numpy won't take it as int32
        if len(self.dds_targets) == 1:
            self._cfr1_ram_configs = [int64(
                (1 << self.bit_shift_ram_enable) |  # ram_enable
                (ad9910.RAM_DEST_ASF << self.bit_shift_ram_dest) |  # ram_destination
                (phase_autoclear << self.bit_shift_phase_autoclear) |
                # (1 << 16) |  # select_sine_output (note: removed for consistency)
                2  # sdio_input_only + msb_first
            ) & 0xFFFFFFFF  # ensure 32b only
            ]
        else:
            self._cfr1_ram_configs += [int64(
                (1 << self.bit_shift_ram_enable) |  # ram_enable
                (ad9910.RAM_DEST_ASF << self.bit_shift_ram_dest) |  # ram_destination
                (phase_autoclear << self.bit_shift_phase_autoclear) |
                # (1 << 16) |  # select_sine_output (note: removed for consistency)
                2  # sdio_input_only + msb_first
            ) & 0xFFFFFFFF # ensure 32b only
                                       ]



    """
    INITIALIZE & CLEAN UP
    """
    @kernel(flags={"fast-math"})
    def sequence_initialize_single_dds(self, dds_targets_idx) -> TNone:
        """
        Prepare experiment hardware for this sequence.
        Should be called when the parent experiment enters the kernel.

        NOTE BENE: this function pulses the CPLD io_update and will affect other dds on the chip
        """
        self.core.break_realtime()
        dds_target = self.dds_targets[dds_targets_idx]

        # disable RAM mode & set matched latencies
        dds_target.set_cfr1(ram_enable=0)
        dds_target.cpld.io_update.pulse_mu(8)
        dds_target.set_cfr2(matched_latency_enable=1)
        dds_target.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

        # write waveform to RAM via RAMWriter
        self.ram_writers[dds_targets_idx].write(self.ampl_asf_pulseshapes_list[dds_targets_idx], self.ram_addr_start)
        delay_mu(25000)

        # set up RAM profile correctly after waveform uploaded
        dds_target.set_profile_ram(
            start=self.ram_addr_start, end=self.ram_addr_stop,
            step=0xFFF,
            profile=self.ram_profile, mode=ad9910.RAM_MODE_RAMPUP
        )
        dds_target.cpld.io_update.pulse_mu(8)
        delay_mu(10000)

        # # clear ASF and POW registers
        dds_target.set_asf(0x00)
        dds_target.set_pow(0x00)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def sequence_initialize(self):
        """
        Initialize all DDSes in the pulse shaper
        """
        for dds_targets_idx in range(len(self.dds_targets)):
            self.sequence_initialize_single_dds(dds_targets_idx)


    @kernel(flags={"fast-math"})
    def sequence_cleanup_single_dds(self, dds_targets_idx) -> TNone:
        """
        Leaves DDS card (specific to this sequence) in a safe state.
        Should be called when the parent experiment completes the main experiment sequence and finishes.
        """
        self.core.break_realtime()
        dds_target = self.dds_targets[dds_targets_idx]

        # stop output & clear registers
        dds_target.sw.off()
        dds_target.set_asf(0x00)
        dds_target.set_pow(0x00)
        dds_target.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

        # disable RAM mode & return CFR1 to normal
        dds_target.set_cfr1(ram_enable=0)
        dds_target.cpld.io_update.pulse_mu(8)
        delay_mu(256)   # add extra slack to avoid RTIO collisions
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def sequence_cleanup(self) -> TNone:
        """
        Leave experiment hardware, i.e all DDSes in the pulse shpaer, (specific to this sequence) in a safe state.
        Should be called when the parent experiment completes the main experiment sequence and finishes.
        """
        for dds_targets_idx in range(len(self.dds_targets)):
            self.sequence_cleanup_single_dds(dds_targets_idx)


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
    def run_single_dds(self, dds_targets_idx,
                       ) -> TNone:
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
        dds_target = self.dds_targets[dds_targets_idx]
        dds_target.set_pow(phas_pow)

        # set RAM profile parameters for amplitude pulse shaping
        dds_target.set_profile_ram(start=self.ram_addr_start, end=self.ram_addr_stop,
                                        profile=self.ram_profile, mode=ad9910.RAM_MODE_RAMPUP,
                                        step=time_ram_step)
        dds_target.cpld.io_update.pulse_mu(8)   # ensure profile is latched

        # set target DDS profile
        dds_target.cpld.set_profile(self.ram_profile)
        dds_target.cpld.io_update.pulse_mu(8)

        # prepare CFRs
        dds_target.write32(ad9910._AD9910_REG_CFR1, int32(self._cfr1_ram_configs[dds_targets_idx])) # enable RAM
        dds_target.cpld.io_update.pulse_mu(8)
        delay_mu(256)   # necessary to prevent RTIO collisions


        ##### FIRE PULSE #####
        time_start_mu = now_mu() & ~7 # coarse align to SYNC_CLK for determinacy
        at_mu(time_start_mu)
        dds_target.cpld.io_update.pulse_mu(8)

        # open and close switch to synchronize with RAM pulse
        at_mu(time_start_mu + self.ram_firing_delay)
        self.external_switches[dds_targets_idx].on()
        dds_target.sw.on()
        delay_mu(time_pulse_mu)
        dds_target.sw.off()
        self.external_switches[dds_targets_idx].off()

        ##### CLEANUP #####
        # disable RAM and DRG
        dds_target.set_cfr1(ram_enable=0)
        dds_target.cpld.io_update.pulse_mu(8)


    """
    TRAIN METHODS
    """
    # todo: note rigor of phase continuity, assumptions made, focus on speed
    @kernel(flags={"fast-math"})
    def configure_train_single_dds(self, dds_targets_idx, time_mu: TInt64) -> TInt64:
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

        dds_target = self.dds_targets[dds_targets_idx]
        # calculate step size/timing for RAM
        time_ram_step = round(time_mu * self.time_pulse_mu_to_ram_step)
        if (time_ram_step > (1 << 16)) or (time_ram_step < 1):
            raise ValueError("Invalid RAM timestep in DDSPulseShape.configure_train().\n"
                             "Change either pulse time or num_samples.")
        # reconvert to get correct time_pulse_mu correctly for later delay
        self.time_pulse_mu_list[dds_targets_idx] = int64(time_ram_step / self.time_pulse_mu_to_ram_step)

        # set RAM profile parameters for pulse shaping
        # note: using RAM rampup mode for simplicity
        dds_target.set_profile_ram(
            start=self.ram_addr_start, end=self.ram_addr_stop,
            step=time_ram_step,
            profile=self.ram_profile, mode=ad9910.RAM_MODE_RAMPUP
        )
        dds_target.cpld.io_update.pulse_mu(8)   # ensure profile is latched

        # switch to target profile
        dds_target.cpld.set_profile(self.ram_profile)
        dds_target.cpld.io_update.pulse_mu(8)   # ensure profile is latched

        # enable RAM mode and leave primed
        dds_target.write32(ad9910._AD9910_REG_CFR1,int32(self._cfr1_ram_configs[dds_targets_idx]))
        dds_target.cpld.io_update.pulse_mu(8)  # ensure profile is latched

        return self.time_pulse_mu_list[dds_targets_idx]

    @kernel(flags={"fast-math"})
    def configure_train_all_dds(self, time_mu_list: TList(TInt64)) -> TList(TInt64):
        """

        """
        if len(time_mu_list) == 1:
            time_mu_list_temp = time_mu_list * len(self.dds_targets)
        elif len(time_mu_list) !=  len(self.dds_targets):
            raise ValueError("Must provide a time array that has the same length as dds_targets.")
        else:
            time_mu_list_temp = time_mu_list

        for dds_targets_idx in range(len(self.dds_targets)):
            self.time_pulse_mu_list[dds_targets_idx] = self.configure_train_single_dds(dds_targets_idx, time_mu_list_temp[dds_targets_idx])

        return self.time_pulse_mu_list


    @kernel(flags={"fast-math"})
    def run_train_single_dds(self, dds_targets_idx: TInt32) -> TNone:
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

        dds_target = self.dds_targets[dds_targets_idx]
        time_start_mu = now_mu() & ~7 # coarse align to SYNC_CLK for determinacy
        at_mu(time_start_mu)
        dds_target.cpld.io_update.pulse_mu(8) # fire pulse!

        # open and close switch to synchronize with RAM pulse
        at_mu(time_start_mu + self.ram_firing_delay)
        dds_target.sw.on()
        delay_mu(self.time_pulse_mu_list[dds_targets_idx])
        dds_target.sw.off()

        # note: we don't clean up (i.e. unset CFR1 or change profiles) b/c we want the
        #   phase accumulator to keep counting - this means that subsequent pulses will
        #   be phase-coherent (i.e. have a deterministic and well-defined phase relationship)
        #   to this pulse.

    @kernel(flags={"fast-math"})
    def run_train_multiple_dds(self, dds_targets_idxs: TList(TInt32)) -> TNone:
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

        time_start_mu = now_mu() & ~7  # coarse align to SYNC_CLK for determinacy
        at_mu(time_start_mu)
        dds_target[dds_targets_idxs[0]].cpld.io_update.pulse_mu(8)  # fire pulse!
        #
        # open and close switch to synchronize with RAM pulse
        at_mu(time_start_mu + self.ram_firing_delay)
        for dds_targets_idx in range(len(self.dds_targets)):
            # open and close switch to synchronize with RAM pulse
            self.external_switches[dds_targets_idx].on()
            self.dds_targets[dds_targets_idx].sw.on()
        for dds_targets_idx in range(len(self.dds_targets)):
            at_mu(time_start_mu + self.ram_firing_delay+ self.time_pulse_mu_list[dds_targets_idx])
            self.dds_targets[dds_targets_idx].sw.off()
            self.external_switches[dds_targets_idx].off()


    @kernel(flags={"fast-math"})
    def run_train_all_dds(self) -> TNone:
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

        time_start_mu = (now_mu() + 8) & ~7  # coarse align to SYNC_CLK for determinacy
        at_mu(time_start_mu)
        self.dds_targets[0].cpld.io_update.pulse_mu(8)  # fire pulse!
        for dds_targets_idx in range(len(self.dds_targets)):
            # open and close switch to synchronize with RAM pulse
            at_mu(time_start_mu + self.ram_firing_delay)
            self.dds_targets[dds_targets_idx].sw.on()
            self.external_switches[dds_targets_idx].on()
        for dds_targets_idx in range(len(self.dds_targets)):
            at_mu(time_start_mu + self.ram_firing_delay + self.time_pulse_mu_list[dds_targets_idx])
            self.dds_targets[dds_targets_idx].sw.off()

    """
    HELPER FUNCTIONS
    """

    class DummyTTL:
        """
        Dummy class so artiq does not get throw an error for a None -type object if no external switch is specified
        """

        @kernel(flags={"fast-math"})
        def on(self):
            delay_mu(8)

        @kernel(flags={"fast-math"})
        def off(self):
            delay_mu(8)
