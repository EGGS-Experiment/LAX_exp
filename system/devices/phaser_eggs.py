from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice
# todo: on/off, set carrier via DUC, pulse shape


# todo: build arguments - att
# todo: prepare - sample time, frame time
# todo: on/off func
# todo: att func


class PhaserEGGS(LAXDevice):
    """
    Device: Phaser (EGGS)

    Use the Phaser AWG to create the EGGS RF
    """
    name = "phaser_eggs"
    core_device = ('phaser', 'phaser0')

    def prepare_device(self):
        # alias devices
        self.phaser_channel =       self.phaser.channel[0]

        # get frequency parameters
        self.freq_cooling_ftw =     self.get_parameter('freq_pump_cooling_mhz', group='beams.freq_mhz', override=False, conversion_function=hz_to_ftw, units=MHz)
        self.freq_readout_ftw =     self.get_parameter('freq_pump_readout_mhz', group='beams.freq_mhz', override=False, conversion_function=hz_to_ftw, units=MHz)
        self.freq_rescue_ftw =      self.get_parameter('freq_pump_rescue_mhz', group='beams.freq_mhz', override=False, conversion_function=hz_to_ftw, units=MHz)

        # get amplitude parameters
        self.ampl_cooling_asf =     self.get_parameter('ampl_pump_cooling_pct', group='beams.ampl_pct', override=False, conversion_function=pct_to_asf)
        self.ampl_readout_asf =     self.get_parameter('ampl_pump_readout_pct', group='beams.ampl_pct', override=False, conversion_function=pct_to_asf)
        self.ampl_rescue_asf =      self.get_parameter('ampl_pump_rescue_pct', group='beams.ampl_pct', override=False, conversion_function=pct_to_asf)

    @kernel(flags={'fast-math'})
    def initialize_device(self):
        # initialize phaser board
        self.phaser.init(debug=True)
        self.core.break_realtime()

        # set nco to center eggs rf around 85 MHz exactly
        self.phaser_channel.set_nco_frequency(-217.083495 * MHz)
        self.phaser_channel.set_nco_phase(0.)
        self.phaser.dac_sync()
        self.core.break_realtime()

        # trf setup, and disable rf output while we set things up
        self.phaser_channel.set_att(0 * dB)
        self.phaser_channel.en_trf_out(rf=0, lo=0)
        self.core.break_realtime()

        # duc
        self.phaser_channel.set_duc_frequency(0 * MHz)
        self.phaser_channel.set_duc_cfg()
        self.phaser.duc_stb()
        self.core.break_realtime()

        # set oscillators (i.e. sidebands)
        self.phaser_channel.oscillator[0].set_frequency_mu(-self.freq_sideband_mu)
        delay_mu(self.time_sample_mu)
        self.phaser_channel.oscillator[0].set_amplitude_phase_mu(0, clr=0)
        delay_mu(self.time_sample_mu)
        self.phaser_channel.oscillator[1].set_frequency_mu(self.freq_sideband_mu)
        delay_mu(self.time_sample_mu)
        self.phaser_channel.oscillator[1].set_amplitude_phase_mu(0, clr=0)
        self.core.break_realtime()

        # re-enable rf output
        self.ph_channel.en_trf_out(rf=1, lo=0)
        self.core.break_realtime()


    @kernel(flags={"fast-math"})
    def set_carrier(self, freq_ftw):
        pass

    @kernel(flags={"fast-math"})
    def set_oscillator_amplitude_phase_mu(self, osc_num, ampl_asf, phase_pow):
        pass


    @kernel(flags={"fast-math"})
    def clear_oscillator(self, osc_num):
        pass
