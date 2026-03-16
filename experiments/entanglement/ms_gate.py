from artiq.experiment import *
from artiq.coredevice import ad9910

from numpy import array, int32, int64, zeros

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, Readout, ReadoutAdaptive,
    RescueIon, QubitRAP
)


class MolmerSorensen(LAXExperiment, Experiment):
    """
    Experiment: Molmer Sorensen Gate


    """
    name = 'Molmer Sorensen Gtae'

    kernel_invariants = {}

    def build_experiment(self):

        # allocate relevant beam profiles
        self.profile_729_RAP = 0
        self.profile_729_SBC = 1
        self.profile_729_readout = 2
        self.profile_729_ms = 3
        self.profile_729_parity = 4

        # devices
        self.setattr_device('qubit')
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')

        self.setattr_argument('reps', NumberValue(default=50, min=10, max=1e6,
                                                  step = 1, scale=, precision=0), tooltip='Number of shots to run')

        self._build_arguements_ms()
        self._build_arguements_parity()


        '''
        Build subsequences
        '''
        self.sidebandcool_subsequence = SidebandCoolContinuousRAM(self, profile_729=self.profile_729_SBC,
        profile_854=3, ram_addr_start_729=0, ram_addr_start_854=0m num_samples=200)
        self.initialize_subsequence = InitializeQubit(self)
        self.readout_subsequence = Readout(self)
        self.rescue_subsequence = RescueIon(self)



    def _build_arguements_ms(self):

        _name = 'ms'

        self.setattr_argument('time_gate_us_list', Scannable(default=[ExplicitScan([50]),
                                                                RangeScan(0.01, 100, 100)],
                                                                global_min=0.01, global_max=1e6,
                                                                global_step = 0.01, unit="us",
                                                                precision=4, scale=1),
                              group= _name, tooltip='time to apply birchromatic interaction for'
                              )


        self.setattr_argument('freq_carrier_ms_mhz',NumberValue(default=101.3,
                                                                min=90, max=101,
                                                                step = 0.01, unit="MHz",
                                                                precision=4, scale=1),
                              group = _name,tooltip='carrier frequency in MHz'
                              )


        self.setattr_argument('detuning_ms_khz_list', Scannable(default=[ExplicitScan([40]),
                                                                RangeScan(0.01, 100, 100)],
                                                                global_min=0.01, global_max=1e6,
                                                                global_step = 0.01, unit="kHz",
                                                                precision=4, scale=1),
                              group = _name,tooltip='detuning from carrier'
                              )

        self.setattr_argument("ampl_ms_pct", PYONValue([50., 50.]), group='default.cat',
                              tooltip="DDS amplitudes for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_ampl_pct, bsb_ampl_pct], which are applied to [singlepass1, singlepass2].")
        self.setattr_argument("atts_ms_db", PYONValue([11., 11.]), group='default.cat',
                              tooltip="DDS attenuations for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_att_db, bsb_att_db], which are applied to [singlepass1, singlepass2].")


    def _build_arguements_parity(self):

        _name = '_parity_pulse'

        self.setattr_argument('enable_parity_pulse', BooleanValue(False),
                              tooltip='Enable parity analysis pulse after applying the MS gate',
                              group=_name)

        self.setattr_argument('time_parity_pulse_us',
                              NumberValue(default=2, min=0.01, max=1e6,
                                          step=0.01, scale=1, ndecimals=4, unit='us'),
                              tooltip='', group=_name)

        self.setattr_argument('phase_parity_pulse_turns_list',
                            Scannable(default=[
                                ExplicitScan([0.]),
                                RangeScan(0, 1, 11),
                                CenterScan(0.5, 1, 0.1)], unit='turns',
                                global_min=0., global_max=2., global_step=0.1, precision=3, scale=1.0),
                            tooltip = '',
                            group = _name)

        self.setattr_argument('att_parity_pulse_db', NumberValue(default=5., min=5., max=31.5,
                                                                step=0.5, ndecimals=1, scale=1, unit='dB'),
                                                                tooltip='', group=_name)

        self.setattr_argument('ampl_parity_pulse_pct', NumberValue(default=50., min=0.01, max=50.,
                                                                step=0.01, ndecimals=3, scale=1, unit='%'),
                                                                tooltip='', group=_name)




    def prepare_experiment(self):
        '''
        Prepare all necessary functions
        '''

        '''
        Convert MS arguments
        '''
        self.time_gate_mu_list = [self.core.seconds_to_mu(time_gate_us*us) for time_gate_us in
                                  self.time_gate_us_list]
        self.freq_carrier_ms_mu = self.qubit.frequency_to_ftw(self.freq_carrier_ms_MHz * MHz)
        self.freq_detuning_ms_mu_list = [self.qubit.frequency_to_ftw(detuning_khz*kHz) for detuning_kHz in
                                    self.freq_detuning_ms_khz_list]
        self.ampl_ms_asf = self.qubit.amplitude_to_asf(self.ampl_ms_pct)


        '''
        Convert Parity Arguments
        '''

        self.time_parity_pulse_us = self.core.seconds_to_mu(self.time_parity_pulse_us*us)
        self.phase_parity_pulse_pow_list = [self.qubit.turns_to_pow(phase_parity_pulse)
                                            for phase_parity_pulse in self.phase_parity_pulse_turns_list]
        self.ampl_parity_pulse_asf = self.qubit.amplitude_to_asf(self.ampl_parity_pulse_pct/100.)



        '''
        Prepare Attenuation Register
        '''

        # attenuation register - bichromatic
        self.att_reg_ms = 0x00000000 | (
                (att_to_mu(self.att_doublepass_default_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (att_to_mu(31.5*dB) << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_ms_db[0] * dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_ms_db[1] * dB) << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

        # attenuation register - parity
        self.att_reg_cat_interferometer = 0x00000000 | (
                (att_to_mu(self.att_doublepass_default_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (att_to_mu(self.att_parity_pulse_db) << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(31.5*dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (att_to_mu(31.5 * dB) << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )



    def initialize_experiment(self) -> TNone:

    def run_main(self) -> TNone:

    def analyze_experiment(self) -> TNone:




