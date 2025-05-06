from numpy import int64
from artiq.experiment import *
from artiq.coredevice.dac34h84 import DAC34H84
from artiq.coredevice.trf372017 import TRF372017

TRF_CONFIG_302_MHZ_CH0 = {
    # LO (i.e. carrier) frequency config
    'pll_div_sel':          0b01,
    'rdiv':                 3,
    'nint':                 29,
    'prsc_sel':             0,
    'cal_clk_sel':          0b1101,

    'icp':                  0b00000,
    'icp_double':           0,

    # freq division config
    'lo_div_sel':           0b11,
    'lo_div_bias':          0b00,
    'bufout_bias':          0b00,

    'tx_div_sel':           0b11,
    'tx_div_bias':          0b00,

    # quadrature modulation compensation (carrier feedthrough)
    # note: for ioff & qoff, mid-range (i.e. 0x88) is zero,
    # 0x00 is most negative, 0xFF is most positive
    'ioff':                 0x00,   # 8b
    'qoff':                 0x84,   # 8b
    'dcoffset_i':           0b01,

    # vco config
    'vco_bias':             0x8,
    'vcobias_rtrim':        0b110,
    'vcobuf_bias':          0b10,

    'vcomux_bias':          0b11,
    'vco_ampl_ctrl':        0b11,
    'vco_vb_ctrl':          0b00,
    'vref_sel':             0b100,

    'pllbias_rtrim':        0b10
}

TRF_CONFIG_302_MHZ_CH1 = {
    # LO (i.e. carrier) frequency config
    'pll_div_sel':          0b01,
    'rdiv':                 3,
    'nint':                 29,
    'prsc_sel':             0,
    'cal_clk_sel':          0b1101,

    'icp':                  0b00000,
    'icp_double':           0,

    # freq division config
    'lo_div_sel':           0b11,
    'lo_div_bias':          0b00,
    'bufout_bias':          0b00,

    'tx_div_sel':           0b11,
    'tx_div_bias':          0b00,

    # quadrature modulation compensation (carrier feedthrough)
    # note: for ioff & qoff, mid-range (i.e. 0x88) is zero,
    # 0x00 is most negative, 0xFF is most positive
    'ioff':                 0xCF,   # 8b
    'qoff':                 0x8F,   # 8b
    'dcoffset_i':           0b10,

    # vco config
    'vco_bias':             0x8,
    'vcobias_rtrim':        0b110,
    'vcobuf_bias':          0b10,

    'vcomux_bias':          0b11,
    'vco_ampl_ctrl':        0b11,
    'vco_vb_ctrl':          0b00,
    'vref_sel':             0b100,

    'pllbias_rtrim':        0b10
}

TRF_CONFIG_781_MHZ_CH0 = {
    # LO (i.e. carrier) frequency config
    'pll_div_sel':          0b01,
    'rdiv':                 2,
    'nint':                 25,
    'prsc_sel':             0,
    'cal_clk_sel':          0b1110,

    'icp':                  0b00000,
    'icp_double':           0,

    # freq division config
    'lo_div_sel':           0b11,
    'lo_div_bias':          0b00,
    'bufout_bias':          0b00,

    'tx_div_sel':           0b10,
    'tx_div_bias':          0b11,

    # quadrature modulation compensation (carrier feedthrough)
    # note: for ioff & qoff, mid-range (i.e. 0x88) is zero,
    # 0x00 is most negative, 0xFF is most positive
    'ioff':                 0x02,   # 8b
    'qoff':                 0x81,   # 8b
    'dcoffset_i':           0b01,

    # vco config
    'vco_bias':             0x8,
    'vcobias_rtrim':        0b110,
    'vcobuf_bias':          0b10,

    'vcomux_bias':          0b11,
    'vco_ampl_ctrl':        0b11,
    'vco_vb_ctrl':          0b00,

    'pllbias_rtrim':        0b10
}

TRF_CONFIG_781_MHZ_CH1 = {
    # LO (i.e. carrier) frequency config
    'pll_div_sel':          0b01,
    'rdiv':                 2,
    'nint':                 25,
    'prsc_sel':             0,
    'cal_clk_sel':          0b1110,

    'icp':                  0b00000,
    'icp_double':           0,

    # freq division config
    'lo_div_sel':           0b11,
    'lo_div_bias':          0b00,
    'bufout_bias':          0b00,

    'tx_div_sel':           0b10,
    'tx_div_bias':          0b11,

    # quadrature modulation compensation (carrier feedthrough)
    # note: for ioff & qoff, mid-range (i.e. 0x88) is zero,
    # 0x00 is most negative, 0xFF is most positive
    'ioff':                 0xCA,   # 8b
    'qoff':                 0x94,   # 8b
    'dcoffset_i':           0b10,

    # vco config
    'vco_bias':             0x8,
    'vcobias_rtrim':        0b110,
    'vcobuf_bias':          0b10,

    'vcomux_bias':          0b11,
    'vco_ampl_ctrl':        0b11,
    'vco_vb_ctrl':          0b00,
    'vref_sel':             0b100,

    'pllbias_rtrim':        0b10
}

DAC_CONFIG_302_MHZ = {
    # general config
    'fifo_offset':          5,
    'invsinc_ena':          0b00,   # MSB is channel AB
    'interpolation':        1,      # x2 interpolation (default)

    # mixer block config
    'mixer_ena':            1,
    'nco_ena':              1,
    'mixer_gain':           1,

    # quadrature modulation compensation (QMC) config
    # QMC gain & phase compensation
    'qmc_corr_ena':         0b11,   # MSB is channel AB

    # note: qmc_gain is [0, 1.9990]; 0x00 = 0, 0x400 = 1.0, 0x7FF = 1.9990
    'qmc_gaina':            0x400,  # 11b (0x7FF max)
    'qmc_gainb':            0x400,  # 11b (0x7FF max)
    'qmc_gainc':            0x40D,  # 11b (0x7FF max)
    'qmc_gaind':            0x3FA,  # 11b (0x7FF max)

    # note: qmc_phase is [-0.5, 0.49975]; two's complement formatted
    'qmc_phaseab':          0x000,  # 12b (default 0x000)
    'qmc_phasecd':          0xFDD,  # 12b (default 0x000)

    # QMC offset compensation
    # note: QMC offset compensation useless here since DAC34H84 is AC-coupled to TRF
    'qmc_offset_ena':       0b00,   # MSB is channel AB
    'qmc_offseta':          0x000,  # 12b
    'qmc_offsetb':          0x000,  # 12b
    'qmc_offsetc':          0x000,  # 12b
    'qmc_offsetd':          0x000,  # 12b
}

DAC_CONFIG_781_MHZ = {
    # general config
    'fifo_offset':          5,
    'invsinc_ena':          0b00,   # MSB is channel AB
    'interpolation':        1,      # x2 interpolation (default)

    # mixer block config
    'mixer_ena':            1,
    'nco_ena':              1,
    'mixer_gain':           1,

    # quadrature modulation compensation (QMC) config
    # QMC gain & phase compensation
    'qmc_corr_ena':         0b00,   # MSB is channel AB

    # note: qmc_gain is [0, 1.9990]; 0x00 = 0, 0x400 = 1.0, 0x7FF = 1.9990
    'qmc_gaina':            0x400,  # 11b (0x7FF max)
    'qmc_gainb':            0x400,  # 11b (0x7FF max)
    # 'qmc_gainc':            0x40D,  # 11b (0x7FF max)
    # 'qmc_gaind':            0x400,  # 11b (0x7FF max)
    'qmc_gainc':            0x3E4,  # 11b (0x7FF max)
    'qmc_gaind':            0x3F4,  # 11b (0x7FF max)

    # note: qmc_phase is [-0.5, 0.49975]; two's complement formatted, max = 0x7FF, min = 0x800, zero = 0x000
    # 'qmc_phaseab':          0x000,  # 12b (default 0x000)
    # 'qmc_phasecd':          0xFDD,  # 12b (default 0x000)
    'qmc_phaseab':          0x000,  # 12b (default 0x000)
    'qmc_phasecd':          0x004,  # 12b (default 0x000)

    # note: grp_delay is [0, t_max]; max = 0xFF, min = zero = 0x00
    'grp_delaya':           0x00,   # 8b (default 0x00)
    'grp_delayb':           0x00,   # 8b (default 0x00)
    'grp_delayc':           0x00,   # 8b (default 0x00)
    'grp_delayd':           0x00,   # 8b (default 0x00)

    # QMC offset compensation
    # note: QMC offset compensation useless here since DAC34H84 is AC-coupled to TRF
    'qmc_offset_ena':       0b00,   # MSB is channel AB
    'qmc_offseta':          0x000,  # 12b
    'qmc_offsetb':          0x000,  # 12b
    'qmc_offsetc':          0x000,  # 12b
    'qmc_offsetd':          0x000,  # 12b
}


class PhaserConfigure(EnvExperiment):
    """
    Tool: Phaser Configure

    Initialize and configure the selected phaser device.
    """

    def build(self):
        # set devices for initialization
        self.setattr_device("core")

        # get list of valid phaser devices and set them as arguments
        phaser_device_list = self._get_phaser_devices()
        self.setattr_argument("phaser_target",      EnumerationValue(list(phaser_device_list), default='phaser1'))

        # frequency configuration
        self.setattr_argument("freq_nco_mhz",       NumberValue(default=0., precision=6, step=100, min=-400., max=400.))
        # self.setattr_argument("freq_trf_mhz",       EnumerationValue(["N/A", "302.083853", "781.251239"], default="302.083853"))
        self.setattr_argument("freq_trf_mhz",       EnumerationValue(["N/A", "302.083853", "781.251239"], default="N/A"))

    def _get_phaser_devices(self):
        """
        Get all valid phaser devices from the device_db.
        """
        def is_local_phaser_device(v):
            return isinstance(v, dict) and (v.get('type') == 'local') and ('class' in v) and (v.get('class') == "Phaser")

        # get only local phaser devices from device_db
        return set([k for k, v in self.get_device_db().items() if is_local_phaser_device(v)])

    def prepare(self):
        """
        Prepare kernel values before running.
        """
        try:
            # get relevant phaser device and add to kernel invariants
            self.phaser = self.get_device(self.phaser_target)
            self.kernel_invariants = self.kernel_invariants | {"phaser"}
        except Exception as e:
            raise Exception("Unable to instantiate target phaser device: {:s}.".format(self.phaser_target))

        # set relevant values for phaser initialization
        self.time_phaser_sample_mu = int64(40)
        self.kernel_invariants = self.kernel_invariants | {"time_phaser_sample_mu"}

        # ensure NCO frequency is valid
        if (self.freq_nco_mhz > 400.) or (self.freq_nco_mhz < -400.):
            raise Exception("Invalid phaser NCO frequency. Must be in range [-400, 400].")
        elif (self.freq_nco_mhz > 300.) or (self.freq_nco_mhz < -300.):
            print("Warning: Phaser NCO frequency outside passband of [-300, 300] MHz.")

        # set up TRF and DAC configuration
        if self.freq_trf_mhz == "N/A":
            self.configure_trf = False
            self.phaser.dac_mmap = DAC34H84().get_mmap()
        else:
            # override phaser object's trf_mmap and dac_mmap for correct operation later on
            self.configure_trf = True
            if self.freq_trf_mhz == "781.251239":
                self.phaser.channel[0].trf_mmap = TRF372017(TRF_CONFIG_781_MHZ_CH0).get_mmap()
                self.phaser.channel[1].trf_mmap = TRF372017(TRF_CONFIG_781_MHZ_CH1).get_mmap()
                self.phaser.dac_mmap = DAC34H84(DAC_CONFIG_781_MHZ).get_mmap()
            elif self.freq_trf_mhz == "302.083853":
                self.phaser.channel[0].trf_mmap = TRF372017(TRF_CONFIG_302_MHZ_CH0).get_mmap()
                self.phaser.channel[1].trf_mmap = TRF372017(TRF_CONFIG_302_MHZ_CH1).get_mmap()
                self.phaser.dac_mmap = DAC34H84(DAC_CONFIG_302_MHZ).get_mmap()
            else:
                raise Exception("Invalid TRF frequency.")

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        """
        Initialize and configure hardware elements on the phaser.
        """
        self.core.reset()

        '''
        *************ATTENUATORS*******************
        '''
        # set maximum attenuation to eliminate output
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_att(31.5 * dB)
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[1].set_att(31.5 * dB)
        delay_mu(10000)


        '''
        *************OSCILLATORS*******************
        '''
        # reset oscillator frequency and amplitude, and keep phase accumulator persistently cleared
        # note: this has to happen before TRF or attenuator adjustment to ensure channel outputs are 0
        for i in range(5):
            # clear channel 0 oscillator
            at_mu(self.phaser.get_next_frame_mu())
            self.phaser.channel[0].oscillator[i].set_frequency(0. * MHz)
            delay_mu(self.time_phaser_sample_mu)
            self.phaser.channel[0].oscillator[i].set_amplitude_phase(amplitude=0., clr=1)

            # clear channel 1 oscillator
            at_mu(self.phaser.get_next_frame_mu())
            self.phaser.channel[1].oscillator[i].set_frequency(0. * MHz)
            delay_mu(self.time_phaser_sample_mu)
            self.phaser.channel[1].oscillator[i].set_amplitude_phase(amplitude=0., clr=1)

        # add slack
        delay_mu(25000)


        '''
        *************PHASER/TRF*******************
        '''
        # configure TRF via phaser initialization (this is the easiest way!)
        if self.configure_trf:
            delay_mu(100000) # 100us
            self.phaser.init(debug=True)

            # ensure TRF outputs are disabled while we finish configurating
            self.core.break_realtime()
            delay_mu(100000) # 100us
            at_mu(self.phaser.get_next_frame_mu())
            self.phaser.channel[0].en_trf_out(rf=0, lo=0)
            delay_mu(self.time_phaser_sample_mu)
            self.phaser.channel[1].en_trf_out(rf=0, lo=0)


        '''
        *************DAC/NCO*******************
        '''
        # set DAC NCO frequency to center output at 85 MHz exactly
        # note: TRF372017 freq is 302.083918 MHz => DAC NCO should be 217.083918 MHz
        # note: currently using -217.083495 as NCO center? idk why
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_nco_frequency(self.freq_nco_mhz * MHz)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser.channel[1].set_nco_frequency(self.freq_nco_mhz * MHz)

        # clear DAC NCO phase
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_nco_phase(0.)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser.channel[1].set_nco_phase(0.)

        # sync DAC for both channels
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.dac_sync()

        # add slack
        self.core.break_realtime()


        '''
        *************DUC (DIGITAL UPCONVERTER)*******************
        '''
        # set channel DUC frequencies to 0
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_duc_frequency(0. * MHz)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser.channel[1].set_duc_frequency(0. * MHz)

        # clear channel DUC phase accumulators
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_duc_cfg(clr_once=1)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser.channel[1].set_duc_cfg(clr_once=1)

        # strobe update register for both DUCs to latch changes
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.duc_stb()

        '''
        *************TRF (UPCONVERTER)*******************
        '''
        # enable outputs for both channels here
        # note: want to do this at end instead of beginning since output may be nonzero and
        #   will be cycling through frequencies as we initialize components
        # note: want to leave trf outputs persistently enabled since phase relation
        #   between channels can change after adjusting the TRF
        if self.configure_trf:
            self.core.break_realtime()
            at_mu(self.phaser.get_next_frame_mu())
            self.phaser.channel[0].en_trf_out(rf=1, lo=0)
            delay_mu(self.time_phaser_sample_mu)
            self.phaser.channel[1].en_trf_out(rf=1, lo=0)

        # add slack
        self.core.break_realtime()


        '''
        **** testing ***
        '''
        # at_mu(self.phaser.get_next_frame_mu())
        # self.phaser.channel[0].set_att_mu(0x00)
        # # self.phaser.channel[0].set_att(0. * dB)
        # at_mu(self.phaser.get_next_frame_mu())
        # # self.phaser.channel[1].set_att_mu(0xFF)
        # self.phaser.channel[1].set_att_mu(0xFF)
        # delay_mu(10000)
        #
        # at_mu(self.phaser.get_next_frame_mu())
        # self.phaser.channel[0].oscillator[0].set_amplitude_phase(amplitude=0.00, clr=0)
        # at_mu(self.phaser.get_next_frame_mu())
        # self.phaser.channel[1].oscillator[0].set_amplitude_phase(amplitude=0.2, clr=0)
        #
        # at_mu(self.phaser.get_next_frame_mu())
        # self.phaser.channel[0].oscillator[0].set_frequency(1. * MHz)
        # at_mu(self.phaser.get_next_frame_mu())
        # self.phaser.channel[1].oscillator[0].set_frequency(1. * MHz)
        # self.core.break_realtime()
        #
        # at_mu(self.phaser.get_next_frame_mu())
        # self.phaser.channel[0].set_duc_frequency(0. * MHz)
        # at_mu(self.phaser.get_next_frame_mu())
        # self.phaser.channel[1].set_duc_frequency(0. * MHz)
        # delay_mu(10000)
        #
        # at_mu(self.phaser.get_next_frame_mu())
        # self.phaser.channel[0].set_duc_cfg(clr=0)
        # at_mu(self.phaser.get_next_frame_mu())
        # self.phaser.channel[1].set_duc_cfg(clr=0)
        # at_mu(self.phaser.get_next_frame_mu())
        # self.phaser.duc_stb()
        # delay_mu(10000)

    def analyze(self):
        """
        Show phaser output center frequency.
        """
        if self.configure_trf:
            freq_trf_mhz = float(self.freq_trf_mhz)
            print("\tPhaser center freq.:\t\t{:.6f} MHz".format(freq_trf_mhz + self.freq_nco_mhz))
            print("\t\tTRF LO Leakage:\t\t{:.6f} MHz".format(freq_trf_mhz))
            print("\t\tTRF SSB Leakage:\t{:.6f} MHz".format(freq_trf_mhz - self.freq_nco_mhz))
            # todo: add spur predict
