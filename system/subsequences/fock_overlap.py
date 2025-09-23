from artiq.experiment import *

from LAX_exp.language import *
from LAX_exp.system.subsequences.qubit_RAP import QubitRAP


class FockOverlap(QubitRAP):
    """
    Subsequence: Fock Overlap

    Inherits from QubitRAP class.
    Use the motional overlap technique from F.Wolf (P.O. Schmidt group) to directly interrogate
        the population of a given fock state (https://www.nature.com/articles/s41467-019-10576-4).
    The maximum interrogable fock state is determined by the number of available auxiliary states.
    Uses RAP for EACH pulse to minimize pulse errors.
    """
    name = 'fock_overlap'
    kernel_invariants = {
        # RAP configuration lists
        "rap_conf_words_list", "rap_conf_att_list", "rap_conf_time_list",

        # RAP sequence specifications
        "sequence_generate", "sequence_readout",
    }

    def build_subsequence(self, ram_profile: TInt32 = -1, ram_addr_start: TInt32 = 0x00,
                          num_samples: TInt32 = 200, pulse_shape: TStr = "blackman"):
        """
        Defines the main interface for the subsequence.
        :param ram_profile: the AD9910 RAM profile to use for RAP + pulse shaping.
        :param ram_addr_start: the beginning RAM register address for pulse shaping.
            Must be in [0, 923].
        :param num_samples: the number of samples to use for the pulse shape.
            Must result in a final RAM address <= 1023.
        :param pulse_shape: the pulse shape to use. Must be supported by available_pulse_shapes.
        """
        # extend our kernel_invariants with parent's kernel_invariants (since we redefine them here)
        kernel_invariants_parent = getattr(super(), "kernel_invariants", set())
        self.kernel_invariants = self.kernel_invariants | kernel_invariants_parent

        _argstr = "fock" # create short string for argument grouping

        # general configuration
        self.setattr_argument("config_rsb",     PYONValue([100.7658, 72., 400., 8.]),
                              group="{}.general".format(_argstr),
                              tooltip="RSB RAP config: [freq_mhz, freq_dev_ss_khz, time_us, att_db].")
        self.setattr_argument("config_bsb",     PYONValue([101.4291, 72., 400., 8.]),
                              group="{}.general".format(_argstr),
                              tooltip="BSB RAP config: [freq_mhz, freq_dev_ss_khz, time_us, att_db].")
        self.setattr_argument("config_carr",    PYONValue([101.0959, 350., 30., 8.]),
                              group="{}.general".format(_argstr),
                              tooltip="Carrier RAP config: [freq_mhz, freq_dev_ss_khz, time_us, att_db].")
        self.setattr_argument("ampl_fock_pct",  NumberValue(default=50., precision=3, step=5., min=0.001, max=50., unit="%", scale=1.),
                              group="{}.general".format(_argstr),
                              tooltip="Same amplitude will be used for all of overlap generation & readout - "
                                      "i.e. RAP, shelving, and carrier pulses.")

        # fock state: generation
        self.setattr_argument("fock_num_prep", NumberValue(default=0, precision=0, step=1, min=0, max=10),
                              group='{}.config'.format(_argstr),
                              tooltip="Fock state number to generate.")

        # fock state: readout
        self.setattr_argument("fock_num_read",  EnumerationValue(["0", "1", "2", "3"], default="0"),
                              group="{}.config".format(_argstr),
                              tooltip="Fock state number to read. Limited by available number of auxiliary states.")
        self.setattr_argument("enable_final_rap", BooleanValue(default=True),
                              group="{}.config".format(_argstr),
                              tooltip="Enable final RAP pulse (which is almost always necessary). "
                                      "This should be ON by default. "
                                      "This option is used for state diagnostics via BSB Rabi divination.")
        self.setattr_argument("config_shelve",  PYONValue({0: (113.2929, 117., 66., 8.),
                                                           1: (104.1451, 60., 274., 8.),
                                                           2: (110.2436, 36., 420., 8.)}),
                              group="{}.config".format(_argstr),
                              tooltip="Shelving RAP config: {fock_num: (freq_mhz, freq_dev_ss_khz, time_us, att_db).")

        # build parent subsequence (QubitRAP)
        super().build_subsequence(
            ram_profile=ram_profile, ram_addr_start=ram_addr_start, num_samples=num_samples,
            ampl_max_pct=self.ampl_fock_pct, pulse_shape=pulse_shape
        )

    def prepare_subsequence(self):
        """
        Prepare and precompute experiment values for speedy evaluation.
        """
        # call parent prepare_subsequence
        super().prepare_subsequence()
        # question: if I do super().prepare_subsequence, does it call their _prepare_argument_checks,
        #  or MY _prepare_argument_checks?
        # answer: it runs MY _prepare_argument_checks (reasonably)


        '''CONVERT RAP CONFIGURATIONS'''
        # need to convert to int since EnumerationValue only allows str
        self.fock_num_read = int(self.fock_num_read)

        # note: handle n=0 case this way b/c idc rubbish params if not fock-ing
        config_shelve_tmp = self.config_shelve if self.fock_num_read != 0 else {0: (85., 43., 20., 8.)} # some dummy config

        # convert all RAP configs
        config_core_all = {'r': self.config_rsb, 'b': self.config_bsb, 'c': self.config_carr, **config_shelve_tmp}
        # note: separate 32b config word list and 32b att from 64b time to avoid conversion & b/c
        #   elements in artiq list must all be of same type
        self.rap_conf_words_list, self.rap_conf_time_list, self.rap_conf_att_list = zip(
            *tuple(
                (*self.configure_values(self.core.seconds_to_mu(conf_dj[2] * us),
                                        self.qubit.frequency_to_ftw(conf_dj[0] * MHz),
                                        self.qubit.frequency_to_ftw(conf_dj[1] * kHz)),
                 att_to_mu(conf_dj[3] * dB))
                for conf_dj in config_core_all.values()
            )
        )

        # note: convert to list b/c artiq doesn't allow tuples to be indexed by variables in-kernel
        self.rap_conf_words_list = list(self.rap_conf_words_list)
        self.rap_conf_time_list = list(self.rap_conf_time_list)
        self.rap_conf_att_list = list(self.rap_conf_att_list)

        # get conversion from sequence keys to index (for later compilation/batch sequencing)
        # config_idx = tuple(range(len(config_core_all)))
        config_idx = dict(zip(tuple(str(i) for i in config_core_all.keys()), range(len(config_core_all))))


        '''COMPILE PULSE SEQUENCES'''
        ### compile generation sequence
        # specify pulses for each n - alternate between RSB and BSB
        seq_generate_str = str.join('r', 'b' * self.fock_num_prep)[: self.fock_num_prep]
        # if odd number of pulses, do a carrier pulse to ensure bring back to S-1/2 to avoid decay event
        seq_generate_str += 'c' if self.fock_num_prep % 2 == 1 else ''
        
        # convert pulse string to parameter index (i.e. location in e.g. rap_conf_words_list)
        if len(seq_generate_str) != 0:
            self.sequence_generate = [config_idx[rap_key] for rap_key in seq_generate_str]
        # account for case where sequence list is empty (artiq forbids empty lists in-kernel)
        else:   self.sequence_generate = [-1]

        ### compile readout sequence
        # specify pulses for each shelving operation (RSB => shelve_<num> => carrier) with final RSB for fluorescence readout
        seq_read_str = ''.join(tuple('r{:d}c'.format(fock_n) for fock_n in range(self.fock_num_read)))
        # configure final RAP according to experiment argument
        seq_read_str += 'r' if self.enable_final_rap else ''
        
        # convert pulse string to parameter index (i.e. location in e.g. rap_conf_words_list)
        if len(seq_read_str) != 0:
            self.sequence_readout = [config_idx[rap_key] for rap_key in seq_read_str]
        # account for case where sequence list is empty (artiq forbids empty lists in-kernel)
        else:   self.sequence_readout = [-1]

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        super()._prepare_argument_checks()

        '''RAP CONFIGS'''
        # check all configs for validity
        config_shelve_tmp = {} if self.fock_num_read == 0 else self.config_shelve
        configs_all = (self.config_rsb, self.config_bsb, self.config_carr, *(config_shelve_tmp.values()))

        for conf_dj in configs_all:
            # check config values are stored correctly
            if not isinstance(conf_dj, (tuple, list)) or (len(conf_dj) != 4):
                raise ValueError("Invalid config: {:}. Must be list or tuple of length 4.".format(conf_dj))
            
            # check config values are of correct type
            elif not all((isinstance(val, (int, float)) for val in conf_dj)):
                raise ValueError("Invalid config values: {:}. Values must all be floats.".format(conf_dj))
            
            # check specific RAP values in range
            elif not (20. <= conf_dj[0] <= 400.): # RAP freq
                raise ValueError("Invalid RAP config freq: {:}. "
                                 "Value must be in range [20, 400] MHz.".format(conf_dj[0]))
            elif conf_dj[1] <= 1.:  # RAP ss freq dev
                raise ValueError("Invalid RAP config freq dev: {:}. "
                                 "Value must be >1 kHz.".format(conf_dj[1]))
            elif conf_dj[2] <= 5.:  # RAP time
                raise ValueError("Invalid RAP config time: {:}. "
                                 "Value must be >5 us.".format(conf_dj[2]))
            elif not (7.9 <= conf_dj[3] <= 31.6): # RAP attenuation
                raise ValueError("Invalid RAP config attenuation: {:}. "
                                 "Value must be in range [8, 31.5] dB.".format(conf_dj[3]))

        '''READOUT'''
        # ensure fock_num_read is valid
        if int(self.fock_num_read) not in (0, 1, 2, 3):
            raise ValueError("Invalid fock_num_read: {:}. "
                             "Must be in ['0', '1', '2', '3'].".format(self.fock_num_read))

        # only require further config_shelve checks if we actually need them for higher |n> detection
        elif self.fock_num_read != "0":
            shelving_keys = self.config_shelve.keys()

            # check config keys valid
            if not all((isinstance(k, int) and (0 <= k <= 2) for k in shelving_keys)):
                raise ValueError("Invalid shelving config keys: {:}. "
                                 "Must be dict with keys as contiguous ints in [0, 2].".format(shelving_keys))
            # ensure shelving config keys are correctly numbered
            elif not ((max(shelving_keys) + 1 == len(shelving_keys)) and (len(set(shelving_keys)) == len(shelving_keys))):
                raise ValueError("Invalid shelving config keys: {:}. "
                                 "Must be contiguous ints starting in [0, 2].".format(shelving_keys))

            # ensure enough shelving for fock state readout
            if max(shelving_keys) + 1 < int(self.fock_num_read):
                raise ValueError("Insufficient shelving configs ({:}) for target fock_num_read ({:}).".format(
                    max(shelving_keys), self.fock_num_read))


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_subsequence(self) -> TNone:
        """
        Prepare the subsequence immediately before run.
        note: copy over parent QubitRAP's intialize_subsequence b/c need to extend it
            and ARTIQ doesn't support in-kernel super().
        """
        '''FROM PARENT (QubitRAP)'''
        # disable RAM + DRG and set matched latencies
        self.qubit.set_cfr1(ram_enable=0)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.qubit.set_cfr2(matched_latency_enable=1)
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

        # write waveform to RAM via RAMWriter
        self.ram_writer.write(self.ampl_asf_pulseshape_list, self.ram_addr_start)
        delay_mu(25000)

    @kernel(flags={"fast-math"})
    def run_fock_generate(self) -> TNone:
        """
        Run fock state generation sequence.
        All pulses are RAP-based.
        For best performance, should be recorded onto DMA.
        """
        self._batch_sequence(self.sequence_generate)

    @kernel(flags={"fast-math"})
    def run_fock_read(self) -> TNone:
        """
        Run fock state generation sequence.
        All pulses are RAP-based.
        For best performance, should be recorded onto DMA.
        """
        self._batch_sequence(self.sequence_readout)

    @kernel(flags={"fast-math"})
    def _batch_sequence(self, sequence_list: TList(TInt32)) -> TNone:
        """
        Run a batch of RAP pulses in succession.
            Used for both fock generation and fock readout b/c all primitives are RAP-based.
        :param sequence_list: list of RAP configurations (denoted by index within the RAP configuration holders,
            e.g. rap_conf_words_list) to be run in sequential order.
        """
        for rap_conf_idx in sequence_list:
            # todo: sanitize input?

            # exit if index value is -1 (workaround b/c we can't have empty lists)
            if rap_conf_idx == -1:
                return

            # set qubit attenuation
            self.qubit.set_att_mu(self.rap_conf_att_list[rap_conf_idx])

            # fire RAP pulse!
            self.run_from_config(self.rap_conf_words_list[rap_conf_idx],
                                 self.rap_conf_time_list[rap_conf_idx])

