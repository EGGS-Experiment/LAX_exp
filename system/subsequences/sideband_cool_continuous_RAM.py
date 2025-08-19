from artiq.experiment import *
from artiq.coredevice import ad9910
from numpy import int32, int64, insert, arange

from LAX_exp.language import *
from LAX_exp.system.objects.RAMWriter import RAMWriter


class SidebandCoolContinuousRAM(LAXSubsequence):
    """
    Subsequence: Sideband Cool - Continuous (RAM)

    Cool the ion to the ground state using a continuous RSB pulse on the S-1/2 to D-5/2 transition.
    Uses RAM mode to reduce the number of profiles taken up and allow greater functionality.
    """
    name = 'sideband_cool_continuous_RAM'
    kernel_invariants = {
        # DDS & SBC config
        "profile_ram_729", "profile_ram_854", "_profile_854_default", "ram_addr_start_729",
        "ram_addr_start_854", "num_samples", "sbc_cycles_cont", "sbc_config_list",

        # beam parameters
        "ampl_qubit_asf", "freq_repump_qubit_ftw", "time_repump_qubit_mu", "time_spinpol_mu", "att_sbc_mu",

        # RAM-related parameters
        "ram_timestep_val", "time_sideband_cooling_mu", "time_spinpolarization_mu_list",
        "ram_waveform_729_ftw_list", "ram_waveform_854_asf_list", "ram_writer_729", "ram_writer_854",
    }

    def build_subsequence(self, profile_729: TInt32 = 1, profile_854: TInt32 = 3,
                          ram_addr_start_729: TInt32 = 0x00, ram_addr_start_854: TInt32 = 0x00,
                          num_samples: TInt32 = 200):
        """
        Defines the main interface for the subsequence.
        Arguments:
            profile_729: the AD9910 RAM profile to use for pulse shaping.
            profile_854: the AD9910 RAM profile to use for pulse shaping.
            ram_addr_start_729: the beginning RAM register address for pulse shaping.
                Must be in [0, 923].
            ram_addr_start_854: the beginning RAM register address for pulse shaping.
                Must be in [0, 923].
            num_samples: the number of samples to use for the pulse shape.
                Must result in a final RAM address <= 1023.
        """
        # set subsequence parameters
        self.profile_ram_729 =      profile_729
        self.profile_ram_854 =      profile_854
        self.ram_addr_start_729 =   ram_addr_start_729
        self.ram_addr_start_854 =   ram_addr_start_854
        self.num_samples =          num_samples

        # get devices
        self.setattr_device('probe')
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('qubit')

        # magic numbers
        self._profile_854_default = 1 # default profile to reset to

        # use helper objects
        self.ram_writer_729 = RAMWriter(self, dds_device=self.qubit.beam,
                                        dds_profile=self.profile_ram_729, block_size=50)
        self.ram_writer_854 = RAMWriter(self, dds_device=self.repump_qubit.beam,
                                        dds_profile=self.profile_ram_854, block_size=50)

    def prepare_subsequence(self):
        """
        Prepare values for speedy evaluation.
        """
        '''PREPARE DATASET PARAMETERS'''
        # sideband cooling configuration parameters
        time_sbc_mu =           self.get_parameter('time_sbc_us', group='sbc', override=False,
                                                   conversion_function=us_to_mu, units=us)
        time_per_spinpol_mu =   self.get_parameter('time_per_spinpol_us', group='sbc', override=False,
                                                   conversion_function=us_to_mu, units=us)
        self.att_sbc_mu =       self.get_parameter('att_sbc_cont_db', group='sbc', override=False,
                                                   conversion_function=att_to_mu, units=dB)

        sbc_cycles_cont =       self.get_parameter('sbc_cycles_cont', group='sbc', override=False)
        self.sbc_config_list =  self.get_parameter('sbc_config_list', group='sbc', override=False)

        # waveform & timing parameters
        self.ampl_qubit_asf =           self.get_parameter('ampl_qubit_pct', group='beams.ampl_pct', override=True,
                                                           conversion_function=pct_to_asf)
        self.freq_repump_qubit_ftw =    self.get_parameter('freq_repump_qubit_mhz', group='beams.freq_mhz', override=False,
                                                           conversion_function=hz_to_ftw, units=MHz)
        self.time_repump_qubit_mu =     self.get_parameter('time_repump_qubit_us', group='timing', override=True,
                                                           conversion_function=seconds_to_mu, units=us)
        self.time_spinpol_mu =          self.get_parameter('time_spinpol_us', group='timing', override=True,
                                                           conversion_function=seconds_to_mu, units=us)

        self._prepare_argument_checks() # note: validate inputs after we get params so they can be checked


        '''PREPARE SIDEBAND COOLING'''
        # calculate steps for each mode, then adjust number of samples to accommodate mode division
        mode_time_pct = tuple(config_arr[0] for config_arr in self.sbc_config_list.values())
        mode_time_steps = tuple(round(val / sum(mode_time_pct) * self.num_samples) for val in mode_time_pct)
        self.num_samples = sum(mode_time_steps)

        # convert timings to multiples of SYNC_CLK (i.e. waveform update clock) period
        time_cycle_mu_to_ram_step = (self.qubit.sysclk_per_mu / 4) / self.num_samples # SYNC_CLK period is 4x AD9910's SYSCLK
        self.ram_timestep_val = round(time_sbc_mu / sbc_cycles_cont * time_cycle_mu_to_ram_step)
        if (self.ram_timestep_val > ((1 << 16) - 1)) or (self.ram_timestep_val < 1):
            raise ValueError("Invalid RAM timestep in SidebandCoolContinuousRAM. Change either num_samples or SBC time.")
        # reconvert to get actual/correct SBC time for later use
        self.time_sideband_cooling_mu = int64(self.ram_timestep_val / time_cycle_mu_to_ram_step * sbc_cycles_cont)

        # calculate spinpol timings
        self.time_spinpolarization_mu_list = arange(
            32, # some small start time number - 32ns
            self.time_sideband_cooling_mu - 100000,   # ensure no spinpol within 100us of SBC end
            time_per_spinpol_mu + self.time_spinpol_mu + 200,  # account for spinpol pulse times & switch delay
            dtype=int64
        )


        '''PREPARE RAM WAVEFORM'''
        # prepare RAMWriters (b/c only LAXExperiment classes call their own children)
        self.ram_writer_729.prepare()
        self.ram_writer_854.prepare()

        # create 729nm waveform array - frequency values
        mode_freqs_hz = tuple(freq_mhz * MHz for freq_mhz in self.sbc_config_list.keys())
        vals_freq_hz = sum(tuple([mode_freqs_hz[i]]*num_steps for i, num_steps in enumerate(mode_time_steps)), [])
        self.ram_waveform_729_ftw_list = [int32(0)] * self.num_samples # stores RAM data in machine units
        self.qubit.frequency_to_ram(vals_freq_hz[::-1], self.ram_waveform_729_ftw_list) # pre-reverse list b/c write_ram reverses it

        # create 854nm waveform array - amplitude values
        mode_quench_ampls_frac = tuple(config_arr[1] / 100. for config_arr in self.sbc_config_list.values())
        vals_ampl_frac = sum(tuple([mode_quench_ampls_frac[i]]*num_steps for i, num_steps in enumerate(mode_time_steps)), [])
        self.ram_waveform_854_asf_list = [int32(0)] * self.num_samples # stores RAM data in machine units
        self.repump_qubit.amplitude_to_ram(vals_ampl_frac[::-1], self.ram_waveform_854_asf_list) # pre-reverse list b/c write_ram reverses it

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # ensure a reasonable amount of modes (i.e. not too many)
        if len(self.sbc_config_list) > 20:
            raise ValueError("Too many modes for SBC. Number of modes must be in [1, 20].")

        # ensure SBC config on all modes add up to 100%
        mode_total_pct = sum(tuple(config_arr[0] for config_arr in self.sbc_config_list.values()))
        if mode_total_pct > 100.:
            raise ValueError("Total sideband cooling mode percentages exceed 100%.")

        # todo: give warning of sum of SBC time % is < 100.
        # todo: ensure total SBC time is in [100, ???] us
        # todo: ensure SBC time step is in [50, ???] us
        # todo: ensure SBC cycle time is in [???, ???] us
        # todo: add checks for these to _prepare_argument_checks
        # time_sbc_mu
        # time_per_spinpol_mu
        # att_sbc_mu
        # sbc_cycles_cont
        # self.sbc_config_list
        # self.ampl_qubit_asf
        # self.freq_repump_qubit_ftw
        # self.time_repump_qubit_mu
        # self.time_spinpol_mu
        # self.setattr_argument("time_sbc_us", NumberValue(default=4000, precision=3, step=100, min=0.001, max=1000000, scale=1., unit='us'))
        # self.setattr_argument("time_per_spinpol_us", NumberValue(default=380, precision=3, step=1, min=0.01, max=100000, scale=1., unit='us'))
        # self.setattr_argument("calib_cont", BooleanValue(default=False))
        # self.setattr_argument("sbc_cycles_cont", NumberValue(default=11, precision=0, step=1, min=1, max=10000))
        # self.setattr_argument("sbc_config_list", PYONValue({100.2832: [37., 6.9], 100.4272: [37., 6.5], 100.7258: [26., 8.0]}))
        # self.setattr_argument("att_sbc_cont_db", NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, scale=1., unit="dB"))


    '''
    COREDEVICE METHODS
    '''
    @kernel(flags={"fast-math"})
    def initialize_subsequence(self) -> TNone:
        """
        Prepare hardware for operation.
        """
        '''729nm RAM SETUP'''
        # disable RAM mode and set matched latencies
        self.qubit.set_cfr1(ram_enable=0)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.qubit.set_cfr2(matched_latency_enable=1)
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

        # write waveform to RAM of 729nm DDS via RAMWriter
        self.ram_writer_729.write(self.ram_waveform_729_ftw_list, self.ram_addr_start_729)
        delay_mu(25000)

        # set up 729nm RAM profile correctly after waveform uploaded
        self.qubit.set_profile_ram(
            start=self.ram_addr_start_729, end=self.ram_addr_start_729 + (self.num_samples - 1),
            step=self.ram_timestep_val,
            profile=self.profile_ram_729, mode=ad9910.RAM_MODE_CONT_RAMPUP
        )
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(10000)


        '''854nm RAM SETUP'''
        # disable RAM mode and set matched latencies
        self.repump_qubit.set_cfr1(ram_enable=0)
        self.repump_qubit.cpld.io_update.pulse_mu(8)
        self.repump_qubit.set_cfr2(matched_latency_enable=1)
        self.repump_qubit.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

        # write waveform to RAM of 729nm DDS via RAMWriter
        self.ram_writer_854.write(self.ram_waveform_854_asf_list, self.ram_addr_start_854)
        delay_mu(25000)

        # set up 854nm RAM profile correctly after waveform uploaded
        self.repump_qubit.set_profile_ram(
            start=self.ram_addr_start_854, end=self.ram_addr_start_854 + (self.num_samples - 1),
            step=self.ram_timestep_val,
            profile=self.profile_ram_854, mode=ad9910.RAM_MODE_CONT_RAMPUP
        )
        self.repump_qubit.cpld.io_update.pulse_mu(8)
        delay_mu(10000)

    @kernel(flags={"fast-math"})
    def cleanup_subsequence(self) -> TNone:
        """
        Clean up the subsequence immediately after run.
        """
        # stop & clear output/registers of SBC beams
        self.qubit.off()
        self.qubit.set_asf(0x00)
        self.qubit.set_ftw(0x00)
        self.qubit.set_pow(0x00)
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

        self.repump_qubit.set_asf(0x00)
        self.repump_qubit.set_ftw(0x00)
        self.repump_qubit.set_pow(0x00)
        self.repump_qubit.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

        # disable RAM mode for SBC beams
        self.qubit.set_cfr1(ram_enable=0)
        self.qubit.cpld.io_update.pulse_mu(8)

        self.repump_qubit.set_cfr1(ram_enable=0)
        self.repump_qubit.cpld.io_update.pulse_mu(8)
        delay_mu(256)   # add extra slack to avoid RTIO collisions

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        """
        Run sideband cooling via RAM mode.
        """
        '''PREPARE HARDWARE FOR SBC '''
        # # tmp remove
        # # quick bugfix to make this subsequence happy (wrt write_ram issues)
        # self.repump_qubit.on()
        # delay_mu(50)
        # self.repump_qubit.off()
        # delay_mu(50)
        # # tmp remove

        with parallel:
            # prepare beams for SBC - 729nm
            with sequential:
                # set SBC waveform
                self.qubit.set_att_mu(self.att_sbc_mu)
                self.qubit.set_asf(self.ampl_qubit_asf)
                self.qubit.cpld.io_update.pulse_mu(8)

                # set target RAM profile & start RAM mode
                self.qubit.cpld.set_profile(self.profile_ram_729)
                self.qubit.cpld.io_update.pulse_mu(8)
                self.qubit.write32(ad9910._AD9910_REG_CFR1,
                                   (1 << 31) |  # ram_enable
                                   (ad9910.RAM_DEST_FTW << 29) |  # ram_destination
                                   (1 << 16) |  # select_sine_output
                                   (1 << 13) |  # phase_autoclear
                                   (1 << 9) | # osk_enable
                                   2 # sdio_input_only + msb_first
                                   )

            # prepare beams for SBC - 854nm
            with sequential:
                # set SBC waveform
                self.repump_qubit.set_ftw(self.freq_repump_qubit_ftw)
                self.repump_qubit.cpld.io_update.pulse_mu(8)

                # set target RAM profile & start RAM mode
                self.repump_qubit.cpld.set_profile(self.profile_ram_854)
                self.repump_qubit.cpld.io_update.pulse_mu(8)
                self.repump_qubit.write32(ad9910._AD9910_REG_CFR1,
                                          (1 << 31) |  # ram_enable
                                          (ad9910.RAM_DEST_ASF << 29) |  # ram_destination
                                          (1 << 16) |  # select_sine_output
                                          (1 << 13) |  # phase_autoclear
                                          2 # sdio_input_only + msb_first
                                          )


        '''BEGIN SBC'''
        time_sbc_init_mu = now_mu()    # get start reference time

        # run spin polarizations according to schedule
        for time_spinpol_mu in self.time_spinpolarization_mu_list:
            at_mu(time_sbc_init_mu + time_spinpol_mu)
            self.spin_polarize()

        # turn on SBC beams (after first spinpol)
        at_mu(time_sbc_init_mu + self.time_spinpol_mu + 10000)
        with parallel:
            with sequential:
                # begin RAM mode and switch on DDS - 729nm
                self.qubit.cpld.io_update.pulse_mu(8)
                self.qubit.on()
            with sequential:
                # begin RAM mode and switch on DDS - 854nm
                self.repump_qubit.cpld.io_update.pulse_mu(8)
                self.repump_qubit.on()

        time_sbc_start_mu = now_mu() # get new fiducial time


        '''FINISH'''
        # stop SBC beams
        at_mu(time_sbc_start_mu + self.time_sideband_cooling_mu)
        with parallel:
            # stop beams via RF switches
            self.qubit.off()
            self.repump_qubit.off()

            # disable RAM mode - 729nm
            with sequential:
                delay_mu(32) # prevents sequence errors
                self.qubit.set_cfr1(ram_enable=0)
                self.qubit.cpld.io_update.pulse_mu(8)

            # disable RAM mode - 854nm
            # note: this is last b/c it takes longest, so timeline is correct after parallel block
            with sequential:
                delay_mu(64) # prevents sequence errors
                self.repump_qubit.set_cfr1(ram_enable=0)
                self.repump_qubit.cpld.io_update.pulse_mu(8)

        # repump qubit after sideband cooling
        delay_mu(32) # prevents sequence errors
        self.repump_qubit.set_profile(self._profile_854_default) # use default (i.e. not SBC quench) profile for 854nm
        self.repump_qubit.cpld.io_update.pulse_mu(8)
        delay_mu(32) # prevents sequence errors
        self.repump_qubit.on()
        delay_mu(self.time_repump_qubit_mu)
        self.repump_qubit.off()

    @kernel(flags={"fast-math"})
    def spin_polarize(self) -> TNone:
        """
        Run spin polarization for optical pumping into the correct Zeeman manifold.
        """
        with parallel:
            self.probe.on()
            self.repump_cooling.on()

        delay_mu(self.time_spinpol_mu)

        with parallel:
            self.probe.off()
            self.repump_cooling.off()

