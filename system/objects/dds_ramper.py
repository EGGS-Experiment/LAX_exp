from artiq.experiment import *
from artiq.coredevice import ad9910, ttl
from numpy import int32, int64, linspace

from LAX_exp.base import LAXDevice

# Digital Ramp Generator - Ramp Destinations
DRG_DEST_FTW =  0b00
DRG_DEST_POW =  0b01
DRG_DEST_ASF =  0b10

class DDSRamper(HasEnvironment):
    """
    Sequence: DDS Ramper

    Apply a coherent, shaped pulse via RAM modulation on the AD9910.
    """
    name = 'DDSRamper'
    kernel_invariants = {
        # objects
        'dds_targets',

        # bit shifters
        '_drg_asf_profile_enable_bit_shift', '_drg_dest_bit_shift', '_drg_enable_bit_shift',
        '_drg_effective_ftw_bit_shift', '_drg_matched_latency_enabled_bit_shift',

        # firing delay
        'ramp_firing_delay',

    }

    def build(self, dds_targets: ad9910.AD9910 = None,
              num_samples: TInt32 = 100,
              ramp_dest: TInt32 = 2,
              data_high: TFloat = 0,
              data_low: TFloat = 0,
              external_switches: ttl.TTLOut = None,
              ramp_time_us: TFloat = 3):
        """
        Interface for using the digital ramp across multiple channels

        :param dds_targets: urukul_channels to use
        :param num_samples: number of samples in the ramp
        :param ramp_dest: ramp destination
            0 ramps frequency, 1 ramps phase and 2 ramps amplitude. All other values will invoke an error
        :param data_high: The high value of the ramp (i.e the maximal value)
        :param data_low: The low value of the ramp (i.e the minimal value)
        :param external_switches: external switches attached to the urukul channel
        :param ramp_time_us: ramp time in microsecond
        """

        # set sequence parameters
        self.drg_ramp_dests = ramp_dest
        self.num_samples = num_samples
        self.external_switches = external_switches

        # get relevant devices
        self.dds_targets = dds_targets
        self.setattr_device("core")

        self.drg_time_ramp_mu = self.core.seconds_to_mu(ramp_time_us*us)

        self._validate_arguments()

        self._cfr1_ram_configs = [0] * len(self.dds_targets)
        self._cfr2_ram_configs = [0] * len(self.dds_targets)

        self.drg_data_high_mu_list = [0] * len(self.dds_targets)
        self.drg_data_low_mu_list = [0] * len(self.dds_targets)
        self.drg_step_mu_list = [0] * len(self.dds_targets)
        self.drg_time_interval_mu_list = [0] * len(self.dds_targets)
        self.drg_time_actual_mu_list = [0] * len(self.dds_targets)

        self._drg_reg_ramp_bit_shift_list = [0] * len(self.dds_targets)

        # bit shifts
        self._drg_asf_profile_enable_bit_shift = 24
        self._drg_dest_bit_shift = 20
        self._drg_enable_bit_shift = 19
        self._drg_effective_ftw_bit_shift = 16
        self._drg_matched_latency_enabled_bit_shift = 7

        # firing delay for ramp
        self.ramp_firing_delay = int64(95)


        ### finish build sequence ###
        for dds_target_idx in range(len(self.dds_targets)):
            self.sequence_prepare(dds_target_idx,
                                  data_high=data_high, data_low=data_low)

    # add auxilary dds target
    def add_dds_target(self,
                       dds_target: ad9910.AD9910 = None,
                       ramp_dest: TInt32 = 2,
                       data_high: TFloat = 0,
                       data_low: TFloat = 0,
                       phase_autoclear: TInt32 = 0,
                       external_switch: ttl.TTLOut = None,
                       ramp_time_us: TFloat = 3) -> TNone:
        """
        Add an additional urukul channel to the digital ramp generator interface

        :param dds_targets: urukul_channels to use
        :param num_samples: number of samples in the ramp
        :param ramp_dest: ramp destination
            0 ramps frequency, 1 ramps phase and 2 ramps amplitude. All other values will invoke an error
        :param data_high: The high value of the ramp (i.e the maximal value)
        :param data_low: The low value of the ramp (i.e the minimal value)
        :param external_switches: external switches attached to the urukul channel
        :param ramp_time_us: ramp time in microsecond
        """

        self.dds_targets += [dds_target]
        self.drg_ramp_dests += [ramp_dest]
        self.drg_time_ramp_mu += [self.core.seconds_to_mu(ramp_time_us*us)]
        self.external_switches += [external_switch]
        self._validate_arguments()

        self.drg_data_high_mu_list += [0]
        self.drg_data_low_mu_list += [0]
        self.drg_step_mu_list += [0]
        self.drg_time_interval_mu_list += [0]
        self.drg_time_actual_mu_list += [0]
        self._drg_reg_ramp_bit_shift_list += [0]


        ### finish build sequence ###
        self.sequence_prepare(len(self.dds_targets) - 1,
                              data_high = data_high, data_low = data_low)

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


        if not isinstance(self.drg_ramp_dests, list):
            self.drg_ramp_dests = [self.drg_ramp_dests]

        if not isinstance(self.drg_time_ramp_mu, list):
            self.drg_time_ramp_mu = [self.drg_time_ramp_mu]


        if not isinstance(self.external_switches, list):
            if self.external_switches is None:
                self.external_switches = [0] * len(self.dds_targets)
            elif isinstance(self.external_switches, ttl.TTLOut):
                self.external_switches = [self.external_switches]

        if not (len(self.dds_targets) == len(
                self.external_switches)):
            raise ValueError("All parameters must have the same shape or be single-valued")

        # ensure dds_target is an AD9910 device
        for dds_target in self.dds_targets:
            if not (isinstance(dds_target, ad9910.AD9910) or isinstance(dds_target, LAXDevice)):
                raise ValueError("Invalid type for dds_target ({}). Must be of type AD9910.".format(
                    type(dds_target)))

        for switch_idx in range(len(self.external_switches)):
            if self.external_switches[switch_idx] == 0 or self.external_switches[switch_idx] is None:
                self.external_switches[switch_idx] = self.DummyTTL()
            elif isinstance(self.external_switches[switch_idx], self.DummyTTL):
                pass
            elif not isinstance(self.external_switches[switch_idx], ttl.TTLOut):
                raise TypeError("the specified the external_switch must be an instance of TTLOut.")

    def sequence_prepare(self, dds_targets_idx: TInt32,
                         data_high: TInt32, data_low: TInt32) -> TNone:
        """
        Prepare & pre-compute relevant values for speedy compilation & kernel evaluation.

        :param dds_targets_idx: index of the urukul channel (within this interface) to prepare
        :param data_high: high value of the ramp (i.e the maximal value)
        :param data_low: low value of the ramp (i.e the minimal value)
        """

        dds_target = self.dds_targets[dds_targets_idx]
        # convert high/low data arguments
        if self.drg_ramp_dests[dds_targets_idx] == DRG_DEST_FTW:
            data_high_mu = dds_target .frequency_to_ftw(data_high)
            data_low_mu = dds_target .frequency_to_ftw(data_low)
            drg_step_mu = dds_target .frequency_to_ftw(
                ((data_high - data_low) / self.num_samples))
        elif self.drg_ramp_dests[dds_targets_idx] == DRG_DEST_POW:
            data_high_mu = dds_target .turns_to_pow(data_high)
            data_low_mu = dds_target .turns_to_pow(data_low)
            drg_step_mu = dds_target .turns_to_pow(
                ((data_high - data_low) / self.num_samples))
        elif self.drg_ramp_dests[dds_targets_idx] == DRG_DEST_ASF:
            data_high_mu = dds_target.amplitude_to_asf(data_high)
            data_low_mu = dds_target.amplitude_to_asf(data_low)
            drg_step_mu = dds_target.amplitude_to_asf(
                ((data_high - data_low) / self.num_samples))
        else:
            raise ValueError("Invalid Ramp destination - must be set to 0 for frequency ramp, 1 for phase ramp, or 2 for "
                             "amplitude ramp")

        if drg_step_mu < 1: #todo check this condition
            drg_step_mu = 1
            # raise ValueError("Need a larger step size. Reduce number of samples or increase ramp range")

        # see ad9910 datasheet
        drg_time_interval =  round(self.drg_time_ramp_mu[dds_targets_idx]*ns * dds_target.sysclk / (4*self.num_samples))
        if drg_time_interval <  1:  # todo: check this condition
            raise ValueError("Need a larger interval time_interval")


        self.drg_data_high_mu_list[dds_targets_idx] = data_high_mu
        self.drg_data_low_mu_list[dds_targets_idx] = data_low_mu
        self.drg_step_mu_list[dds_targets_idx] = drg_step_mu
        self.drg_time_interval_mu_list[dds_targets_idx] = round(drg_time_interval)
        self.drg_time_actual_mu_list[dds_targets_idx] = int64(drg_time_interval * self.num_samples)

        if self.drg_ramp_dests[dds_targets_idx] == DRG_DEST_FTW:
            self._drg_reg_ramp_bit_shift_list[dds_targets_idx] = 0
        elif self.drg_ramp_dests[dds_targets_idx] == DRG_DEST_POW:
            self._drg_reg_ramp_bit_shift_list[dds_targets_idx] = 16
        else:
            self._drg_reg_ramp_bit_shift_list[dds_targets_idx] = 18


    """
    INITIALIZE & CLEAN UP
    """

    @kernel(flags={"fast-math"})
    def sequence_initialize_single_dds(self, dds_targets_idx: TInt32) -> TNone:
        """
        Prepare experiment hardware for this sequence.
        Should be called when the parent experiment enters the kernel.
        """
        self.core.break_realtime()
        dds_target = self.dds_targets[dds_targets_idx]

        # set basic configuration
        dds_target.set_cfr1()
        dds_target.set_cfr2()


    @kernel(flags={"fast-math"})
    def sequence_initialize(self):
        """
        Initialize all DDSes in the ramp generator interface
        """
        for dds_targets_idx in range(len(self.dds_targets)):
            self.sequence_initialize_single_dds(dds_targets_idx)

    @kernel(flags={"fast-math"})
    def sequence_cleanup_single_dds(self, dds_targets_idx) -> TNone:
        """
        Leaves DDS card (specific to this sequence) in a safe state.
        Should be called when the parent experiment completes the main experiment sequence and finishes.

        :param dds_targets_idx: index of the urukul channel (within this interface)

        NOTE BENE: this function pulse io_update and will potentially effect other channels on the chip
        """
        self.core.break_realtime()
        dds_target = self.dds_targets[dds_targets_idx]

        # stop output & clear registers
        dds_target.sw.off()
        delay_mu(25000)

        # disable RAM mode & return CFR1 and CFR2 to normal
        dds_target.set_cfr1()
        dds_target.set_cfr2(matched_latency_enable=1)
        dds_target.cpld.io_update.pulse_mu(8)
        delay_mu(256)  # add extra slack to avoid RTIO collisions
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
    TRAIN METHODS
    """

    # todo: note rigor of phase continuity, assumptions made, focus on speed
    @kernel(flags={"fast-math"})
    def configure_ramp_single_dds(self, dds_targets_idx: TInt32,
                                  ramp_direction: TInt32 = 1) -> TNone:
        """
        Configure the cfr and ramp registers with precomputed values to enable ramp generation on the next io_update

        :param dds_targets_idx: index of the urukul channel (within this interface)
        :param ramp_direction: rising or falling mode
            0 indicates a falling linear ramp and 1 indicates a rising linear ramp
=
        """
        dds_target = self.dds_targets[dds_targets_idx]

        # configure CFRs
        CFR1_DRG_CONFIG = int64(
            (0 << 13) |  # phase_autoclear
            (1 << 15) |  # load_lrr (this bit pulsed high mean every io_update will start a ramp from our current value
            2  # sdio_input_only + msb_first
        ) & 0xFFFFFFFF  # ensure 32b only
        dds_target.write32(ad9910._AD9910_REG_CFR1, int32(CFR1_DRG_CONFIG)) # enable digital ramp generation
        CFR2_DRG_CONFIG = int64(
            (1 << self._drg_asf_profile_enable_bit_shift) | # asf_profile_enable
            (1 << self._drg_effective_ftw_bit_shift) | # effective_ftw
            (1 << self._drg_matched_latency_enabled_bit_shift) |  # matched_latency_enable
            (self.drg_ramp_dests[dds_targets_idx] << 20) |  # digital_ramp_destination
            (1 << self._drg_enable_bit_shift)  # digital_ramp_enable
        ) & 0xFFFFFFFF # ensure 32b only
        dds_target.write32(ad9910._AD9910_REG_CFR2, int32(CFR2_DRG_CONFIG)) # enable digital ramp generation


        reg_ramp_limit_bit_shift = self._drg_reg_ramp_bit_shift_list[dds_targets_idx]
        dds_target.write64(ad9910._AD9910_REG_RAMP_LIMIT,
                           data_high=self.drg_data_high_mu_list[dds_targets_idx] << reg_ramp_limit_bit_shift,
                           data_low=self.drg_data_low_mu_list[dds_targets_idx] << reg_ramp_limit_bit_shift)

        # set Digital Ramp Generator update interval
        dds_target.write32(ad9910._AD9910_REG_RAMP_RATE,
                           (self.drg_time_interval_mu_list[dds_targets_idx] << 16) |   # ramp down
                           (self.drg_time_interval_mu_list[dds_targets_idx] << 0)    # ramp up
                           )

        if ramp_direction == 1:
            self.configure_rising_ramp_single_dds(dds_targets_idx)
        else:
            self.configure_falling_ramp_single_dds(dds_targets_idx)


    @kernel(flags={"fast-math"})
    def configure_rising_ramp_single_dds(self, dds_targets_idx: TInt32) -> TNone:
        """
        Configure a rising linear ramp

        :param dds_targets_idx: index of the urukul channel (within this interface)
        """
        dds_target = self.dds_targets[dds_targets_idx]
        reg_ramp_limit_bit_shift = self._drg_reg_ramp_bit_shift_list[dds_targets_idx]

        dds_target.write64(ad9910._AD9910_REG_RAMP_STEP,
                           data_high=-self.drg_step_mu_list[dds_targets_idx] << reg_ramp_limit_bit_shift,  # ramp down
                           data_low=self.drg_step_mu_list[dds_targets_idx] << reg_ramp_limit_bit_shift  # ramp up
                           )

    @kernel(flags={"fast-math"})
    def configure_rising_ramp_all_dds(self) -> TNone:
        """
        Configure a rising linear ramp for all DDSes in interface
        """
        for dds_targets_idx in range(len(self.dds_targets)):
            self.configure_rising_ramp_single_dds(dds_targets_idx)

    @kernel(flags={"fast-math"})
    def configure_falling_ramp_single_dds(self, dds_targets_idx: TInt32) -> TNone:
        """
        Configure a falling linear ramp

        :param dds_targets_idx: index of the urukul channel (within this interface)
        """
        dds_target = self.dds_targets[dds_targets_idx] 
        reg_ramp_limit_bit_shift = self._drg_reg_ramp_bit_shift_list[dds_targets_idx]

        dds_target.write64(ad9910._AD9910_REG_RAMP_STEP,
                           data_high=self.drg_step_mu_list[dds_targets_idx] << reg_ramp_limit_bit_shift,
                           data_low=-self.drg_step_mu_list[dds_targets_idx] << reg_ramp_limit_bit_shift,
                           )

    @kernel(flags={"fast-math"})
    def configure_falling_ramp_all_dds(self) -> TNone:
        """
        Configure a rising linear ramp for all DDSes in interface
        """
        for dds_targets_idx in range(len(self.dds_targets)):
            self.configure_falling_ramp_single_dds(dds_targets_idx)


    @kernel(flags={"fast-math"})
    def configure_ramp_all_dds(self,
                               ramp_direction) -> TNone:
        """
        Configures the linear ramp for all DDSes in this interface

        :param ramp_direction: direction of ramp
        """
        for dds_targets_idx in range(len(self.dds_targets)):
            self.configure_ramp_single_dds(dds_targets_idx, ramp_direction)


    @kernel(flags={"fast-math"})
    def run_ramp_all_dds(self) -> TNone:
        """
        Fire linear ramped pulse

        This only has an effect if we are not at the maximal value (i.e. running a rising ramp while we are at the high
        value has no effect and we just stay on the high value) OR drg autoclear is set high
        """
        time_start_mu = (now_mu() + 8) & ~7  # coarse align to SYNC_CLK for determinacy
        at_mu(time_start_mu)
        self.dds_targets[0].cpld.io_update.pulse_mu(8)  # fire pulse!
        delay_mu(self.ramp_firing_delay)
        for dds_targets_idx in range(len(self.dds_targets)):
            # open and close switch to synchronize with RAM pulse
            at_mu(time_start_mu)
            self.dds_targets[dds_targets_idx].sw.on()
            self.external_switches[dds_targets_idx].on()

    '''
    Helper functions
    '''
    @kernel(flags={"fast-math"})
    def reset_cfr1_all_dds(self):
        """
        Reset cfr1 register to normal values for all DDSes in this interface
        """
        for dds_target in self.dds_targets:
            dds_target.set_cfr1()

    @kernel(flags={"fast-math"})
    def reset_cfr2_all_dds(self):
        """
        Reset cfr2 register to normal values for all DDSes in this interface
        """
        for dds_target in self.dds_targets:
            dds_target.set_cfr2()

    @kernel(flags={"fast-math"})
    def reset_cfrs_all_dds(self):
        """
        Reset cfr1 and cfr2 register to normal values for all DDSes in this interface
        """
        for dds_target in self.dds_targets:
            dds_target.set_cfr1()
            dds_target.set_cfr2()

    @kernel(flags={"fast-math"})
    def switch_on_all_dds(self):
        """
        Toggle on all DDSes in this interface
        """
        for dds_target_idx in range(len(self.dds_targets)):
            self.dds_targets[dds_target_idx].sw.on()
            delay_mu(8)

    @kernel(flags={"fast-math"})
    def switch_off_all_dds(self):
        """
        Toggle off all DDSes in this interface
        """
        for dds_target_idx in range(len(self.dds_targets)):
            self.dds_targets[dds_target_idx].sw.off()
            delay_mu(8)

    @kernel(flags={"fast-math"})
    def set_autoclear_phase_accumulator_all_dds(self):
        """
        Clear the dds phase accumulator for all DDSes in this interface
        """
        for dds_target in self.dds_targets:
            dds_target.set_cfr1(phase_autoclear=1)

    """
    HELPER CLASSES
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
