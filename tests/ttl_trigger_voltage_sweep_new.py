import labrad
import numpy as np

from os import environ
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE
from EGGS_labrad.config.dc_config import dc_config


class TTLTriggerVoltageSweepNew(EnvExperiment):
    """
    TTL Trigger Voltage SweepNew
    """
    kernel_invariants = {
        'time_timeout_pmt_mu',
        'dc_micromotion_channel',
        'ampl_mod_vpp',
        'freq_mod_mhz',
        'dc_micromotion_voltages_v_list'
    }

    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge"
    ]


    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # repetitions
        self.setattr_argument("repetitions",                        NumberValue(default=20000, ndecimals=0, step=1, min=1, max=10000000))

        # timing
        self.setattr_argument("time_timeout_pmt_ms",                NumberValue(default=10000, ndecimals=5, step=1, min=1, max=1000000))

        # modulation
        self.setattr_argument("ampl_mod_vpp",                       NumberValue(default=0.5, ndecimals=3, step=1, min=1, max=1000000))
        self.setattr_argument("freq_mod_mhz",                       NumberValue(default=1.59, ndecimals=5, step=1, min=0, max=1000000))

        # voltage values
        self.dc_micromotion_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_micromotion_channel",             EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='V Shim'))
        self.setattr_argument("dc_micromotion_voltages_v_list",     Scannable(
                                                                        default=RangeScan(40.0, 80.0, 41, randomize=True),
                                                                        global_min=0, global_max=1000, global_step=1,
                                                                        unit="V", scale=1, ndecimals=4
                                                                    ))

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)

        # connect to labrad
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.fg = self.cxn.function_generator_server
        self.dc = self.cxn.dc_server


    def prepare(self):
        # PMT devices
        self.pmt_counter =                                          self.get_device("ttl{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                      getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_timeout_pmt_mu =                                  self.core.seconds_to_mu(self.time_timeout_pmt_ms * ms)

        # get voltage parameters
        self.dc_micromotion_voltages_v_list =                       np.array(list(self.dc_micromotion_voltages_v_list))
        self.dc_micromotion_channel =                               self.dc_micromotion_channeldict[self.dc_micromotion_channel]['num']

        # RF modulation synchronization clock
        self.mod_clock =                                            self.get_device("urukul0_ch3")
        self.mod_toggle =                                           self.get_device("ttl8")
        self.mod_clock_freq_ftw =                                   self.mod_clock.frequency_to_ftw(10. * MHz)
        self.mod_clock_ampl_pct =                                   self.mod_clock.amplitude_to_asf(0.5)
        self.mod_clock_att_db =                                     2. * dB
        self.mod_clock_delay_mu =                                   self.core.seconds_to_mu(300 * ns)

        # set up datasets
        self._dataset_counter                                       = 0
        self.set_dataset("ttl_trigger",                             np.zeros([len(self.dc_micromotion_voltages_v_list), self.repetitions]))
        self.setattr_dataset("ttl_trigger")

        # record parameters
        self.set_dataset('xArr',                                    self.dc_micromotion_voltages_v_list)
        self.set_dataset('repetitions',                             self.repetitions)
        self.set_dataset('freq_mod_mhz',                            self.freq_mod_mhz)
        self.set_dataset('modulation_amplitude_vpp',                self.ampl_mod_vpp)
        self.set_dataset('dc_channel_num',                          self.dc_micromotion_channel)

        # set up function generator
        # get list of function generators
        fg_dev_list = self.fg.list_devices()
        fg_dev_dict = dict(tuple(fg_dev_list))

        # select correct function generator
        dev_exists = False
        for dev_num, dev_desc in fg_dev_dict.items():
            if 'DG2P' in dev_desc:
                dev_exists = True
                self.fg.select_device(dev_num)

        # raise error if function generator doesn't exist
        if not dev_exists:
            raise Exception("Error: modulation function generator not detected.")

        # set up function generator
        self.fg.toggle(0)
        self.fg.amplitude(self.ampl_mod_vpp)
        self.fg.frequency(self.freq_mod_mhz * 1e6)

        # set up gated burst
        self.fg.burst(True)
        self.fg.burst_mode('GAT')

        # enable output
        self.fg.toggle(1)


    @kernel(flags={"fast-math"})
    def run(self):
        self.core.reset()

        # set ttl directions
        self.pmt_counter.input()
        self.mod_toggle.output()
        self.mod_toggle.off()
        self.core.break_realtime()

        # configure rf mod clock
        self.mod_clock.set_phase_mode(PHASE_MODE_ABSOLUTE)
        self.core.break_realtime()

        # start rf mod clock
        self.mod_clock.set_att(self.mod_clock_att_db)
        self.mod_clock.set_mu(self.mod_clock_freq_ftw, asf=self.mod_clock_ampl_pct)
        self.mod_clock.cfg_sw(True)
        self.core.break_realtime()

        # enable external clocking
        self.fg_write(':ROSC:SOUR EXT')

        # MAIN LOOP
        # sweep voltage
        for voltage_val in self.dc_micromotion_voltages_v_list:

            # reset FIFOs
            self.core.reset()

            # set frequency
            self.voltage_set(self.dc_micromotion_channel, voltage_val)
            self.core.break_realtime()

            # set up loop variables
            counter = 0
            time_stop_mu_list = [0] * self.repetitions
            self.core.break_realtime()

            # synchronize timings with DDS clock
            #self.mod_clock.set_mu(self.mod_clock_freq_ftw, asf=self.mod_clock_ampl_pct)
            #delay_mu(self.mod_clock_delay_mu)

            # activate modulation and wait for change in output
            self.mod_toggle.on()
            delay_mu(self.mod_clock_delay_mu)
            time_start_mu = now_mu()

            # start counting photons
            time_stop_mu = self.pmt_counter.gate_rising_mu(self.time_timeout_pmt_mu)
            while counter < self.repetitions:

                # move timestamped photons
                time_mu_tmp = self.pmt_counter.timestamp_mu(time_stop_mu)

                # increase gating time
                if time_mu_tmp < 0:
                    print('\tError: increased gating time')
                    self.core.break_realtime()
                    time_stop_mu = self.pmt_counter.gate_rising_mu(self.time_timeout_pmt_mu)
                else:
                    time_stop_mu_list[counter] = time_mu_tmp
                    counter += 1

            # sync delay at end
            self.core.break_realtime()

            # complete a loop
            self.mod_toggle.off()
            self.update_dataset(time_start_mu, time_stop_mu_list)
            self.core.break_realtime()


    # LABRAD FUNCTIONS
    @rpc
    def voltage_set(self, channel, voltage_v):
        """
        Set the channel to the desired voltage.
        """
        voltage_set_v = self.dc.voltage(channel, voltage_v)
        print('\tvoltage set: {}'.format(voltage_set_v))

    @rpc
    def frequency_set(self, freq_hz):
        """
        Set the RF to the desired frequency.
        """
        freq_set_hz = self.fg.frequency(freq_hz)
        print('\tfrequency set: {}'.format(freq_set_hz))

    def fg_write(self, msg):
        """
        write a GPIB message todo: clean up
        """
        self.fg.gpib_write(msg)

    @rpc(flags={"async"})
    def update_dataset(self, time_start_mu, time_stop_mu_list):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.mutate_dataset(
            'ttl_trigger',
            self._dataset_counter,
            np.array(self.core.mu_to_seconds(np.array(time_stop_mu_list) - time_start_mu))
        )

        self._dataset_counter += 1

    def analyze(self):
        # turn off modulation
        self.fg.toggle(0)

        # process data
        # ttl_trigger_tmp = np.array(self.ttl_trigger).reshape((len(self.dc_micromotion_voltages_v_list), self.repetitions, 2))
        # ind_arr = np.argsort(self.freq_mod_mhz_list)
        # ttl_trigger_tmp = ttl_trigger_tmp[ind_arr]
        #
        # self.ttl_trigger_processed = ttl_trigger_tmp
