from numpy import int64
from artiq.experiment import *
from artiq.coredevice.trf372017 import TRF372017


class PhaserConfigure(EnvExperiment):
    """
    Utility: Phaser Configure

    Initialize and configure the selected phaser device.
    """

    def build(self):
        # set devices for initialization
        self.setattr_device("core")

        # get list of valid phaser devices and set them as arguments
        phaser_device_list = self._get_phaser_devices()
        self.setattr_argument("phaser_target",      EnumerationValue(list(phaser_device_list), default='phaser1'))

        # frequency configuration
        self.setattr_argument("freq_nco_mhz",       NumberValue(default=-217.083495, ndecimals=6, step=100, min=-400., max=400.))
        self.setattr_argument("freq_trf_mhz",       EnumerationValue(["N/A", 302.083918, 781.251239], default="N/A"))

        # dataset management
        # todo: dataset updating - boolean: freq_center
        # self.setattr_argument("calibration",        BooleanValue(default=False))

    def _get_phaser_devices(self):
        """
        Get all valid phaser devices from the device_db.
        """
        def is_local_phaser_device(v):
            return isinstance(v, dict) and (v.get('type') == 'local') and ('class' in v) and (v.get('class') == "Phaser")

        # get only local phaser devices from device_db
        # todo: make it return values as well so we can get the TRF dict
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
            print("Error: unable to instantiate target phaser device.")
            raise e

        # set relevant values for phaser initialization
        self.time_phaser_sample_mu = int64(40)
        self.kernel_invariants = self.kernel_invariants | {"time_phaser_sample_mu"}

        # ensure NCO frequency is valid
        if (self.freq_nco_mhz > 400.) or (self.freq_nco_mhz < -400.):
            raise Exception("Invalid phaser NCO frequency. Must be in range [-400, 400].")
        elif (self.freq_nco_mhz > 300.) or (self.freq_nco_mhz < -300.):
            print("Warning: Phaser NCO frequency outside passband of [-300, 300] MHz.")

        # set up TRF configuration
        if self.freq_trf_mhz == "None":
            self.configure_trf = False
            self.configure_trf_mmap = []
        else:
            self.configure_trf = True
            if self.freq_trf_mhz == 781.251239:
                trf_config_update = {
                    'rdiv': 2,
                    'nint': 25,
                    'pll_div_sel': 0b01,
                    'prsc_sel': 0,

                    'icp': 0b00000,
                    'icp_double': 0,

                    'cal_clk_sel': 0b1110,

                    'lo_div_sel': 0b11,
                    'lo_div_bias': 0b00,
                    'bufout_bias': 0b00,

                    'tx_div_sel': 0b10,
                    'tx_div_bias': 0b00
                }
            elif self.freq_trf_mhz == 302.083918:
                trf_config_update = {
                    'pll_div_sel':          0b01,
                    'rdiv':                 3,
                    'nint':                 29,
                    'prsc_sel':             0,
                    'cal_clk_sel':          0b1101,

                    'lo_div_sel':           0b11,
                    'lo_div_bias':          0b00,
                    'bufout_bias':          0b00,

                    'tx_div_sel':           0b11,
                    'tx_div_bias':          0b00
                }
            else:
                raise Exception("Invalid TRF frequency.")

            # create a TRF object and get mmap
            trf_object = TRF372017(trf_config_update)
            self.configure_trf_mmap = trf_object.get_mmap()

    @kernel(flags={"fast-math"})
    def run(self):
        """
        todo:document
        :return:
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

        # add slack
        self.core.break_realtime()


        '''
        *************PHASER*******************
        '''
        # initialize phaser
        delay_mu(1000000)
        self.phaser.init(debug=True)

        # add slack
        self.core.break_realtime()

        # ensure TRF outputs are disabled while we configure things
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

        # add slack
        self.core.break_realtime()


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
        self.core.break_realtime()


        '''
        *************TRF (UPCONVERTER)*******************
        '''
        # update TRFs on both channels with mmap to set output freq center
        if self.configure_trf is True:
            delay(1000000)
            for i in range(2):
                for data in self.configure_trf_mmap:
                    self.core.break_realtime()
                    at_mu(self.phaser.get_next_frame_mu())
                    self.phaser.channel[i].trf_write(data)

        # enable outputs for both channels here
        # note: want to do this at end instead of beginning since output may be nonzero and
        #   will be cycling through frequencies as we initialize components
        # note: want to leave trf outputs persistently enabled since phase relation
        #   between channels can change after adjusting the TRF
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].en_trf_out(rf=1, lo=0)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser.channel[1].en_trf_out(rf=1, lo=0)

        # add slack
        self.core.break_realtime()

        # ensure completion
        self.core.wait_until_mu(now_mu())

