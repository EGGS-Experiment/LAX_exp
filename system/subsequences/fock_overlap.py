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
        "config_idx",
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
        _argstr = "fock" # create short string for argument grouping
        self._config_holder_gen = [] # store user-compiled sequence configs for state generation
        self._config_holder_read = [] # store user-compiled sequence configs for state readout

        # extend our kernel_invariants with parent's kernel_invariants (since we redefine them here)
        kernel_invariants_parent = getattr(super(), "kernel_invariants", set())
        self.kernel_invariants = self.kernel_invariants | kernel_invariants_parent

        self.setattr_argument("ampl_fock_pct",  NumberValue(default=50., precision=3, step=5., min=0.001, max=50., unit="%", scale=1.),
                              group="{}.config".format(_argstr),
                              tooltip="Same amplitude will be used for all of overlap generation & readout.\n"
                                      "i.e. for RAP, carrier, and shelving pulses.")
        self.setattr_argument("enable_final_rap", BooleanValue(default=True),
                              group="{}.config".format(_argstr),
                              tooltip="Enable final RAP pulse (which is almost always necessary). "
                                      "This should be ON by default. "
                                      "This option is used for state diagnostics via BSB Rabi divination.")

        # build parent subsequence (QubitRAP)
        super().build_subsequence(
            ram_profile=ram_profile, ram_addr_start=ram_addr_start, num_samples=num_samples,
            ampl_max_pct=self.ampl_fock_pct, pulse_shape=pulse_shape
        )

    def prepare_subsequence(self):
        """
        Prepare and precompute experiment values for speedy evaluation.
        """
        '''
        GENERAL SETUP
        '''
        # get relevant configs from dataset manager
        # note: has to happen before super().prepare_subsequence since (trailed off - need to figure out what I was going to write)
        self.setattr_parameter('config_shelve', group='sequences.fock_overlap')
        self.setattr_parameter('config_rsb', group='sequences.fock_overlap')
        self.setattr_parameter('config_bsb', group='sequences.fock_overlap')
        self.setattr_parameter('config_carrier', group='sequences.fock_overlap')

        # call parent prepare_subsequence
        # note: parent prepare_subsequence runs MY _prepare_argument_checks (reasonably)
        super().prepare_subsequence()


        '''
        CONVERT RAP CONFIGURATIONS
        '''
        config_core_all = {'r': self.config_rsb, 'b': self.config_bsb, 'c': self.config_carrier, **self.config_shelve}
        # note: separate 32b word list and 32b att from 64b time b/c artiq list elements must all be of same type
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
        self.config_idx = dict(
            zip(
                tuple(str(i) for i in config_core_all.keys()),
                range(len(config_core_all))
            )
        )

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        super()._prepare_argument_checks()

        '''RAP CONFIGS'''
        # check all configs for validity
        config_shelve_tmp = self.config_shelve
        configs_all = (self.config_rsb, self.config_bsb, self.config_carrier, *(config_shelve_tmp.values()))

        for conf_dj in configs_all:
            # check config values are stored correctly
            if not isinstance(conf_dj, (tuple, list)) or (len(conf_dj) != 4):
                raise ValueError("Invalid config: {:}. Must be list or tuple of length 4.".format(conf_dj))

            # check config values are of correct type
            elif not all((isinstance(val, (int, float)) for val in conf_dj)):
                raise ValueError("Invalid config values: {:}. Values must all be floats.".format(conf_dj))

            # check specific RAP values in range
            elif not (20. <= conf_dj[0] <= 400.):  # RAP freq
                raise ValueError("Invalid RAP config freq: {:}. "
                                 "Value must be in range [20, 400] MHz.".format(conf_dj[0]))
            elif conf_dj[1] <= 1.:  # RAP ss freq dev
                raise ValueError("Invalid RAP config freq dev: {:}. "
                                 "Value must be >1 kHz.".format(conf_dj[1]))
            elif conf_dj[2] <= 5.:  # RAP time
                raise ValueError("Invalid RAP config time: {:}. "
                                 "Value must be >5 us.".format(conf_dj[2]))
            elif not (7.9 <= conf_dj[3] <= 31.6):  # RAP attenuation
                raise ValueError("Invalid RAP config attenuation: {:}. "
                                 "Value must be in range [8, 31.5] dB.".format(conf_dj[3]))

        '''READOUT'''
        # only require further config_shelve checks if we actually need them for higher |n> detection
        shelving_keys = self.config_shelve.keys()

        # check config keys valid
        if not all((isinstance(k, int) and (0 <= k <= 2) for k in shelving_keys)):
            raise ValueError("Invalid shelving config keys: {:}. "
                             "Must be dict with keys as contiguous ints in [0, 2].".format(shelving_keys))
        # ensure shelving config keys are correctly numbered
        elif not ((max(shelving_keys) + 1 == len(shelving_keys)) and (len(set(shelving_keys)) == len(shelving_keys))):
            raise ValueError("Invalid shelving config keys: {:}. "
                             "Must be contiguous ints starting in [0, 2].".format(shelving_keys))

        # todo: reimplement some checks
        # # ensure enough shelving for fock state readout
        # if max(shelving_keys) + 1 < int(self.fock_num_read):
        #     raise ValueError("Insufficient shelving configs ({:}) for target fock_num_read ({:}).".format(
        #         max(shelving_keys), self.fock_num_read))


    '''
    USER INTERFACE
    '''
    @rpc
    def create_generate_sequence(self, fock_num: TInt32) -> TList(TInt32):
        """
        Prepare compilation of a fock state generation sequence.
        :param fock_num: the fock state number to generate.
        """
        # sanitize input: int from [0, N]
        if not (0 <= fock_num <= 10):
            raise ValueError("Invalid fock_num for create_generate_sequence. Must be int in [0, 10].")

        ### compile generation sequence
        # specify pulses for each n - alternate between RSB and BSB
        seq_generate_str = str.join('r', 'b' * fock_num)[: fock_num]
        # if odd number of pulses, do a carrier pulse to ensure bring back to S-1/2 to avoid decay event
        seq_generate_str += 'c' if fock_num % 2 == 1 else ''

        # convert pulse string to parameter index (i.e. location in e.g. rap_conf_words_list)
        if len(seq_generate_str) != 0:
            return [self.config_idx[rap_key] for rap_key in seq_generate_str]
        # account for case where sequence list is empty (artiq forbids empty lists in-kernel)
        else:   return [-1]

    @rpc
    def create_read_sequence(self, fock_num: TInt32) -> TList(TInt32):
        """
        Prepare compilation of a fock overlap readout sequence.
        :param fock_num: the fock state number to generate.
        """
        # sanitize input: int from [0, 3]
        if not (0 <= fock_num <= 3):
            raise ValueError("Invalid fock_num for create_read_sequence. Must be int in [0, 3].")

        ### compile readout sequence
        # specify pulses for each shelving operation (RSB => shelve_<num> => carrier) with final RSB for fluorescence readout
        seq_read_str = ''.join(tuple('r{:d}c'.format(fock_n) for fock_n in range(fock_num)))
        # configure final RAP according to experiment argument
        seq_read_str += 'r' if self.enable_final_rap else ''

        # convert pulse string to parameter index (i.e. location in e.g. rap_conf_words_list)
        if len(seq_read_str) != 0:
            return [self.config_idx[rap_key] for rap_key in seq_read_str]
        # account for case where sequence list is empty (artiq forbids empty lists in-kernel)
        else: return [-1]


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_subsequence(self) -> TNone:
        """
        Prepare the subsequence immediately before run.
        note: copy over parent QubitRAP's initialize_subsequence b/c need to extend it
            and ARTIQ doesn't support in-kernel super().
        """
        ### FROM PARENT (QubitRAP) ###
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
    def run_sequence(self, sequence_list: TList(TInt32)) -> TNone:
        """
        Run a batch of RAP pulses in succession.
            Used for both fock generation and fock readout b/c all primitives are RAP-based.
        :param sequence_list: list of RAP configurations (denoted by index within the RAP configuration holders,
            e.g. rap_conf_words_list) to be run in sequential order.
        """
        for rap_conf_idx in sequence_list:
            # todo: sanitize input?
            # exit if index value is -1 (workaround b/c we can't have empty lists)
            if rap_conf_idx == -1: return

            # set qubit attenuation
            self.qubit.set_att_mu(self.rap_conf_att_list[rap_conf_idx])

            # fire RAP pulse!
            self.run_from_config(self.rap_conf_words_list[rap_conf_idx],
                                 self.rap_conf_time_list[rap_conf_idx])

