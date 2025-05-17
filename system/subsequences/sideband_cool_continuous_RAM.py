from artiq.experiment import *
from artiq.coredevice import ad9910

import numpy as np
from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class SidebandCoolContinuousRAM(LAXSubsequence):
    """
    Subsequence: Sideband Cool - Continuous (RAM)

    Cool the ion to the ground state using a continuous RSB pulse on the S-1/2 to D-5/2 transition.
    Uses RAM mode to reduce the number of profiles taken up and allow greater functionality.
    """
    name = 'sideband_cool_continuous_RAM'
    kernel_invariants = {
        # general subsequence parameters
        "profile_ram_729", "profile_ram_854", "ram_addr_start_729", "ram_addr_start_854", "num_samples",
        "ampl_qubit_asf", "freq_repump_qubit_ftw", "time_repump_qubit_mu", "time_spinpol_mu",
        "att_sidebandcooling_mu",

        # RAM-related parameters
        "freq_dds_sync_clk_hz", "time_cycle_mu_to_ram_step", "ram_timestep_val", "time_sideband_cooling_mu",
        "time_spinpolarization_mu_list",
        "ram_waveform_729_ftw_list", "ram_waveform_854_asf_list"
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

        # sideband cooling - hardware values
        self.setattr_argument("time_sideband_cooling_us",   NumberValue(default=5000, precision=3, step=100, min=0.001, max=1000000), group='SBC_RAM.continuous')
        self.setattr_argument("time_per_spinpol_us",        NumberValue(default=600, precision=3, step=1, min=0.01, max=100000), group='SBC_RAM.continuous',
                                                                tooltip="time between spin polarization pulses (in us)")

        # sideband cooling - configuration
        self.setattr_argument("calibration_continuous",     BooleanValue(default=False), group='SBC_RAM.continuous',
                                                                tooltip="True: disables 729nm DDS during SBC for calibration purposes")
        self.setattr_argument("sideband_cycles_continuous", NumberValue(default=10, precision=0, step=1, min=1, max=10000), group='SBC_RAM.continuous',
                                                                tooltip="number of times to loop over the SBC configuration sequence")
        self.setattr_argument("sideband_cooling_config_list",       PYONValue({100.7555: [26., 5.], 100.455: [37., 5.], 100.315: [37., 5.]}),
                              group='SBC_RAM.continuous',
                              tooltip="{freq_mode_mhz: [sbc_mode_pct_per_cycle, ampl_quench_mode_pct]}")
        self.setattr_argument("att_sidebandcooling_continuous_db",  NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5), group='SBC_RAM.continuous')

    def prepare_subsequence(self):
        """
        Prepare values for speedy evaluation.
        """
        '''VALIDATE INPUTS'''
        self._prepare_argument_checks()

        '''RETRIEVE PARAMETERS'''
        # waveform parameters
        self.ampl_qubit_asf =       self.get_parameter('ampl_qubit_pct', group='beams.ampl_pct', override=True,
                                                       conversion_function=pct_to_asf)
        self.freq_repump_qubit_ftw = self.get_parameter('freq_repump_qubit_mhz', group='beams.freq_mhz', override=False,
                                                        conversion_function=hz_to_ftw, units=MHz)
        # timing parameters
        self.time_repump_qubit_mu = self.get_parameter('time_repump_qubit_us', group='timing', override=True,
                                                       conversion_function=seconds_to_mu, units=us)
        self.time_spinpol_mu =      self.get_parameter('time_spinpol_us', group='timing', override=True,
                                                       conversion_function=seconds_to_mu, units=us)

        '''CONVERT ARGUMENTS TO MACHINE UNITS'''
        # waveform values
        self.att_sidebandcooling_mu =   att_to_mu(self.att_sidebandcooling_continuous_db * dB)

        # SBC configuration arrays
        mode_freqs_hz = np.array([
            freq_mhz * MHz
            for freq_mhz in self.sideband_cooling_config_list.keys()
        ])
        mode_time_pct = np.array([
            config_arr[0]
            for config_arr in self.sideband_cooling_config_list.values()
        ])
        mode_quench_ampls_frac = np.array([
            config_arr[1] / 100.
            for config_arr in self.sideband_cooling_config_list.values()
        ])

        '''PREPARE SIDEBAND COOLING'''
        # calculate steps for each mode, then adjust number of samples to accommodate mode division
        mode_time_steps = np.int32([
            round(val)
            for val in mode_time_pct / np.sum(mode_time_pct) * self.num_samples
        ])
        self.num_samples = np.sum(mode_time_steps)

        # convert pulse times into multiples of the SYNC_CLK (i.e. waveform update clock) period
        # todo: get sync_clk from ad9910 device instead
        self.freq_dds_sync_clk_hz = 1e9 / 4.  # SYNC_CLK HAS 4ns PERIOD
        self.time_cycle_mu_to_ram_step = ((self.freq_dds_sync_clk_hz / self.core.seconds_to_mu(1)) /
                                          self.num_samples)
        # calculate DDS register value to set the timestep
        self.ram_timestep_val = round(self.core.seconds_to_mu(self.time_sideband_cooling_us * us) /
                                 self.sideband_cycles_continuous * self.time_cycle_mu_to_ram_step)
        if (self.ram_timestep_val > ((1 << 16) - 1)) or (self.ram_timestep_val < 1):
            raise ValueError("Invalid RAM timestep in SidebandCoolContinuousRAM."
                             "Change either number of samples or adjust SBC time.")
        # reconvert to get actual/correct SBC time for later use
        self.time_sideband_cooling_mu = np.int64(self.ram_timestep_val / self.time_cycle_mu_to_ram_step *
                                                 self.sideband_cycles_continuous)

        # calculate spinpol timings
        time_per_spinpol_mu = self.core.seconds_to_mu(self.time_per_spinpol_us * us)
        # note: we account for the nonzero spinpol time in these delays
        self.time_spinpolarization_mu_list = np.arange(
            time_per_spinpol_mu,
            self.time_sideband_cooling_mu - 10000,   # ensure we don't spinpol within 10us of SBC end
            time_per_spinpol_mu + self.time_spinpol_mu + 200,  # add 200ns to account for switches
            dtype=np.int64
        )
        # ensure we do a spinpol at least once at the beginning
        # self.time_spinpolarization_mu_list = np.insert(self.time_spinpolarization_mu_list, 0, 16)

        '''PREPARE RAM WAVEFORM'''
        # create 729nm waveform array - frequency values
        vals_freq_hz = np.concatenate([
            mode_freqs_hz[i] * np.ones(num_steps)
            for i, num_steps in enumerate(mode_time_steps)
        ])
        # create holder to store RAM waveform data in machine units (ftw, 729nm SBC)
        self.ram_waveform_729_ftw_list = [np.int32(0)] * self.num_samples
        # pre-reverse waveform list since write_ram makes a booboo and reverses the array
        self.qubit.frequency_to_ram(vals_freq_hz[::-1], self.ram_waveform_729_ftw_list)

        # create 854nm waveform array - amplitude values
        vals_ampl_frac = np.concatenate([
            mode_quench_ampls_frac[i] * np.ones(num_steps)
            for i, num_steps in enumerate(mode_time_steps)
        ])
        # create holder to store RAM waveform data in machine units (asf, 854nm quench)
        self.ram_waveform_854_asf_list = [np.int32(0)] * self.num_samples
        # pre-reverse waveform list since write_ram makes a booboo and reverses the array
        self.repump_qubit.amplitude_to_ram(vals_ampl_frac[::-1], self.ram_waveform_854_asf_list)

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # ensure a reasonable amount of modes (i.e. not too many)
        if len(self.sideband_cooling_config_list) > 20:
            raise ValueError("Too many modes for SBC. Number of modes must be in [1, 20].")

        # ensure SBC config on all modes add up to 100%
        mode_total_pct = np.sum([config_arr[0] for config_arr in self.sideband_cooling_config_list.values()])
        if mode_total_pct > 100.:
            raise ValueError("Total sideband cooling mode percentages exceed 100%.")

        # todo: ensure total SBC time is in [100, ???] us
        # todo: ensure SBC time step is in [50, ???] us
        # todo: ensure SBC cycle time is in [???, ???] us

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self) -> TNone:
        """
        Prepare hardware for operation.
        """
        self.core.break_realtime()

        # disable RAM mode and set matched latencies
        with parallel:
            with sequential:
                self.qubit.set_cfr1(ram_enable=0)
                self.qubit.set_cfr2(matched_latency_enable=1)
                self.qubit.cpld.io_update.pulse_mu(8)

            with sequential:
                self.repump_qubit.set_cfr1(ram_enable=0)
                self.repump_qubit.set_cfr2(matched_latency_enable=1)
                self.repump_qubit.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # configure RAM waveform profiles - 729nm for SBC
        # prepare to write waveform to RAM profile
        self.qubit.set_profile_ram(
            start=self.ram_addr_start_729, end=self.ram_addr_start_729 + (self.num_samples - 1),
            step=self.ram_timestep_val,
            profile=self.profile_ram_729, mode=ad9910.RAM_MODE_CONT_RAMPUP
        )

        # set target RAM profile
        self.qubit.cpld.set_profile(self.profile_ram_729)
        self.qubit.cpld.io_update.pulse_mu(8)

        # write waveform to RAM profile
        self.core.break_realtime()
        delay_mu(10000000)   # 10 ms
        self.qubit.write_ram(self.ram_waveform_729_ftw_list)
        self.core.break_realtime()

        # configure RAM waveform profiles - 854nm for quench
        # prepare to write waveform to RAM profile
        self.repump_qubit.set_profile_ram(
            start=self.ram_addr_start_854, end=self.ram_addr_start_854 + (self.num_samples - 1),
            step=self.ram_timestep_val,
            profile=self.profile_ram_854, mode=ad9910.RAM_MODE_CONT_RAMPUP
        )

        # set target RAM profile
        self.repump_qubit.cpld.set_profile(self.profile_ram_854)
        self.repump_qubit.cpld.io_update.pulse_mu(8)

        # write waveform to RAM profile
        self.core.break_realtime()
        delay_mu(10000000)  # 10 ms
        self.repump_qubit.write_ram(self.ram_waveform_854_asf_list)
        self.core.break_realtime()
        delay_mu(1000000)

        # bugfix: needed to make RAM/DMA/whatever happy
        self.repump_qubit.on()
        delay_mu(1000)
        self.repump_qubit.off()

    @kernel(flags={"fast-math"})
    def cleanup_subsequence(self) -> TNone:
        """
        Clean up the subsequence immediately after run.
        """
        self.core.break_realtime()

        # stop & clear output/registers of SBC beams
        self.qubit.off()
        self.qubit.set_asf(0x00)
        self.qubit.set_ftw(0x00)
        self.qubit.set_pow(0x00)
        self.qubit.cpld.io_update.pulse_mu(8)

        self.repump_qubit.set_asf(0x00)
        self.repump_qubit.set_ftw(0x00)
        self.repump_qubit.set_pow(0x00)
        self.repump_qubit.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # disable RAM mode for SBC beams
        self.qubit.set_cfr1(ram_enable=0)
        self.qubit.cpld.io_update.pulse_mu(8)

        self.repump_qubit.set_cfr1(ram_enable=0)
        self.repump_qubit.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # add extra slack following cleanup
        delay_mu(100000)   # 100 us
        self.repump_qubit.cpld.io_update.pulse_mu(8)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        """
        Run sideband cooling via RAM mode.
        """
        '''PREPARE ION FOR SBC '''
        # quick bugfix to make this subsequence happy (wrt write_ram issues)
        self.repump_qubit.on()
        delay_mu(50)
        self.repump_qubit.off()
        delay_mu(50)

        # prepare beams for sideband cooling
        with parallel:
            # set sideband cooling attenuation for qubit beam
            with sequential:
                self.qubit.set_att_mu(self.att_sidebandcooling_mu)
                self.qubit.set_asf(self.ampl_qubit_asf)
                self.qubit.cpld.io_update.pulse_mu(8)

            # set sideband cooling profiles for regular beams
            with sequential:
                self.repump_qubit.set_profile(self.profile_ram_854)
                self.repump_qubit.set_ftw(self.freq_repump_qubit_ftw)
                # do spin polarization before SBC (per Guggemos' thesis)
                # self.spin_polarize()

        # '''SCHEDULE SPINPOL'''
        # # note: we do this here due to difficulties w/empty list for spinpol scheduling
        # # get start reference time
        # time_start_mu = now_mu()
        #
        # # do spin polarizations according to schedule
        # for time_spinpol_mu in self.time_spinpolarization_mu_list:
        #     at_mu(time_start_mu + time_spinpol_mu)
        #     self.spin_polarize()

        '''PRIME RAM MODE FOR SBC BEAMS'''
        # at_mu(time_start_mu + self.time_spinpol_mu + 10000)
        with parallel:
            with sequential:
                # set target RAM profile
                self.qubit.cpld.set_profile(self.profile_ram_729)
                self.qubit.cpld.io_update.pulse_mu(8)
                # start RAM mode
                self.qubit.write32(ad9910._AD9910_REG_CFR1,
                                   (1 << 31) |  # ram_enable
                                   (ad9910.RAM_DEST_FTW << 29) |  # ram_destination
                                   (1 << 16) |  # select_sine_output
                                   (1 << 13) |  # phase_autoclear
                                   (1 << 9) | # osk_enable
                                   2
                                   )
                self.qubit.cpld.io_update.pulse_mu(8)   # note: this starts stepping through RAM

                # FIRE BEAM
                if self.calibration_continuous:
                    self.qubit.off()
                else:
                    self.qubit.on()

            with sequential:
                # set target RAM profile
                self.repump_qubit.cpld.set_profile(self.profile_ram_854)
                self.repump_qubit.cpld.io_update.pulse_mu(8)
                # start RAM mode
                self.repump_qubit.write32(ad9910._AD9910_REG_CFR1,
                                          (1 << 31) |  # ram_enable
                                          (ad9910.RAM_DEST_ASF << 29) |  # ram_destination
                                          (1 << 16) |  # select_sine_output
                                          (1 << 13) |  # phase_autoclear
                                          2
                                          )
                self.repump_qubit.cpld.io_update.pulse_mu(8)    # note: this starts stepping through RAM
                self.repump_qubit.on()  # FIRE BEAM


        '''SCHEDULE SPINPOL'''
        # note: we do this here due to difficulties w/empty list for spinpol scheduling
        # get start reference time
        time_start_mu = now_mu()

        # do spin polarizations according to schedule
        for time_spinpol_mu in self.time_spinpolarization_mu_list:
            at_mu(time_start_mu + time_spinpol_mu)
            self.spin_polarize()

        '''FINISH'''
        # stop sideband cooling
        # at_mu(time_start_mu + self.time_spinpol_mu + 200 + self.time_sideband_cooling_mu)
        at_mu(time_start_mu + self.time_sideband_cooling_mu)
        with parallel:
            # stop beams via RF switches
            self.qubit.off()
            self.repump_qubit.off()

            # disable RAM mode - 729nm
            with sequential:
                self.qubit.set_cfr1(ram_enable=0)
                self.qubit.cpld.io_update.pulse_mu(8)

            # disable RAM mode - 854nm
            with sequential:
                self.repump_qubit.set_cfr1(ram_enable=0)
                self.repump_qubit.cpld.io_update.pulse_mu(8)

        # repump qubit after sideband cooling
        # use normal beam profile for 854 (i.e. not SBC quench)
        self.repump_qubit.set_profile(1)
        self.repump_qubit.cpld.io_update.pulse_mu(8)
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

