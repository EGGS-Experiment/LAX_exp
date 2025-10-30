from artiq.experiment import *
from numpy import mean, zeros, array, std

import labrad
from os import environ
from time import time, sleep
from EGGS_labrad.clients import createTrunk


class LabjackRead(EnvExperiment):
    """
    Tool: Labjack Read

    Read values over time from the labjack T-7 ADC.
    """
    name = "LabjackRead"
    kernel_invariants = {
        # measurement configuration
        'channel_list', 'range_list', 'resolution_list',
        'time_start_s', '_dataset_idx',

        # labrad cxn
        'cxn', 'dv', 'labjack', 'c_record',
    }

    def build(self):
        # channels & gains
        self.setattr_argument("channel_range_resolution_dict",  PYONValue({"AIN0": [10, 12], "AIN1": [0.01, 12]}),
                              tooltip="A dictionary of {channel_name: [measurement_range_volts, resolution_bits]}.\n"
                                      "Channel name must be one of [AIN0, ..., AIN13].\n"
                                      "Measurement range is must be a number of [10, 1, 0.1, 0.01].\n"
                                      "Resolution (in bits) must be an int in [0, 12].\n"
                                      "Note that the actual measurement range is bipoler, e.g. a value of 10 sets a range of [-10, +10]V\n"
                                      "Can record multiple channels, though this will limit the max sample rate (due to timing overheads).")

        # timing
        self.setattr_argument("poll_interval_s",    NumberValue(default=0.5, precision=2, step=1, min=0.2, max=100, scale=1., unit="s"),
                              tooltip="Amount of time to wait between subsequent samples.\n")
        self.setattr_argument("num_samples",        NumberValue(default=100, precision=0, step=10, min=1, max=10000000, scale=1., unit="samples"),
                              tooltip="The total number of samples to record for.")

        # set default scheduling so we can run in background
        self.setattr_device("scheduler")
        self.set_default_scheduling(pipeline_name="labrad")

    def prepare(self) -> TNone:
        """
        Prepare relevant values for experiment runtime.
        """
        ### PROCESS INPUTS ###
        # parse input channel configuration
        self.channel_list = [
            port_name.strip().upper()
            for port_name in self.channel_range_resolution_dict.keys()
        ]
        self.range_list, self.resolution_list = [
            list(channel_param_vals)
            for channel_param_vals in zip(*[
                channel_params
                for channel_params in self.channel_range_resolution_dict.values()
            ])
        ]

        # convert channel configs to numerical values
        self.range_list = [float(val) for val in self.range_list]
        self.resolution_list = [int(val) for val in self.resolution_list]

        # validate inputs
        self._prepare_argument_checks()

        ### CREATE DATASETS ###
        self.set_dataset('results', zeros([self.num_samples, 1 + len(self.channel_list)], dtype=float))
        self.setattr_dataset('results')
        self._dataset_idx = 0

    def _prepare_argument_checks(self):
        """
        Check experiment arguments for validity.
        """
        # check channel names
        channel_names_valid = tuple(
            isinstance(port_name, str) and ("AIN" in port_name) and (port_name.split('AIN')[1].isnumeric())
            for port_name in self.channel_list
        )
        if not all(channel_names_valid):
            raise ValueError("Invalid channel name.\n"
                             "All channels must be str 'AIN<port_num>' where port_num is an int in [0, 14].")

        # check channel configurations
        channel_gains_valid = tuple(
            isinstance(gain_val, (int, float)) and
            any(abs(gain_val - gain_target) <= 1e-5 for gain_target in (10, 1, 0.1, 0.01))
            for gain_val in self.range_list
        )
        if not all(channel_gains_valid):
            raise ValueError("Invalid channel gain. Must be one of [10, 1, 0.1, 0.01].")

        channel_res_valid = tuple(
            isinstance(res_val, int) and (1 <= res_val <= 12)
            for res_val in self.resolution_list
        )
        if not all(channel_res_valid):
            raise ValueError("Invalid channel gain. Must be one of [10, 1, 0.1, 0.01].")


    """
    MAIN SEQUENCE
    """
    @rpc
    def _run_prepare(self) -> TNone:
        """
        Prepare relevant hardware only immediately before starting the experiment
            to prevent interfering with other experiments during e.g. prepare.
        """
        ### ESTABLISH LABRAD CONNECTIONS ###
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.labjack = self.cxn.labjack_server
        self.dv = self.cxn.data_vault
        self.c_record = self.cxn.context()

        ### SET UP DATAVAULT FOR REAL-TIME VIEWING (VIA GRAPHER) ###
        self.time_start_s = time() # get start time
        self.dv.cd(createTrunk(self.name), True, context=self.c_record)
        self.dv.new(
            'Labjack Read', # dataset title
            [('Elapsed time', 't')],    # x-axis/independent variable
            [   # y-axis/dependent variables
                ('{:s}'.format(port_name), 'Volts', 'V')
                for port_name in self.channel_list
            ],
            context=self.c_record
        )

        ### SET UP HARDWARE ###
        # configure range & resolution for labjack analog channels
        self.labjack.write_names(   # configure input range
            ["{:s}_RANGE".format(port_name) for port_name in self.channel_list], # address names
            # ['FLOAT32'] * len(self.channel_list), # register data types
            self.range_list # target values
        )
        self.labjack.write_names(   # configure resolution (in bits)
            ["{:s}_RESOLUTION_INDEX".format(port_name) for port_name in self.channel_list], # address names
            # ['UINT16'] * len(self.channel_list), # register data types
            self.resolution_list # target values
        )

    @rpc
    def run(self) -> TNone:
        """
        Run main sampling loop.
        """
        try:
            self._run_prepare() # prepare relevant devices etc

            for i in range(self.num_samples):
                # sample and record data
                self.update_results(i, [self.labjack.read_name(port_name) for port_name in self.channel_list])

                # periodically check termination
                if (i % 5 == 0) and self.scheduler.check_termination(): break

                # wait until next sample
                sleep(self.poll_interval_s)

        finally:
            self.cxn.disconnect()   # clean up

    @rpc
    def update_results(self, i: TInt32, res_arr: TList(TFloat)) -> TNone:
        """
        Save results to relevant datasets.
        :param i: the sample number
        :param res_arr: list of recorded voltages
        """
        time_elapsed_s = time() - self.time_start_s
        self.mutate_dataset("results", i, array([time_elapsed_s] + res_arr))   # save to artiq dataset
        self.dv.add(time_elapsed_s, *res_arr, context=self.c_record)    # save to labrad dataset
        self._dataset_idx += 1  # increment counter


    '''
    ANALYSIS
    '''
    def analyze(self) -> TNone:
        """
        Print summary statistics of results.
        """
        print('\tResults:')
        for idx_chan, channel_name in enumerate(self.channel_list):
            print('\t\t{:s}:\t{:.4g} +/- {:.4g} mV'.format(
                self.channel_list[idx_chan],
                mean(self.results[:self._dataset_idx, 1+idx_chan]) * 1e3,
                std(self.results[:self._dataset_idx, 1+idx_chan]) * 1e3)
            )

