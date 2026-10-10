import numpy as np
from artiq.experiment import *
from matplotlib.pyplot import imsave

import os
from time import time, sleep
from datetime import datetime

from LAX_exp.language import *

#todo: add status checks for F70 to report any errors

class QuickWarmup(LAXExperiment, Experiment):
    """
    Utility: Quick Warmup
    """
    name = 'Quick Warmup'

    def build_experiment(self):
        # general arguments
        group_name_timing = 'timing'
        self.setattr_argument('time_warmup_min',
                              NumberValue(default=4., min=0.5,
                                          max=10., step=0.1, precision=1,
                                          scale=1, unit='min'),
                              group = group_name_timing
                              )

        self.setattr_argument('time_record_after_warmup_min',
                              NumberValue(default=15., min=0.5,
                                          max=60., step=0.1, precision=1,
                                          scale=1, unit='min'),
                              group = group_name_timing
                              )


        # get necessary devices
        self.setattr_device('f70_compressor')
        self.setattr_device('ion_pump')
        self.setattr_device('lakeshore')


    @property
    def results_shape(self):
        return (self.repetitions_quick_warmup + self.repetitions_record_after_warmup,
                6)

    def prepare_experiment(self):
        """
        Prepare experimental values and precompute/preallocate to reduce kernel overheads.
        """
        self.time_warmup_s = self.time_warmup_min * 60
        self.iterer_quick_warmup = np.arange(1, self.time_warmup_s, 1)
        self.repetitions_quick_warmup = len(self.iterer_quick_warmup)

        self.time_record_after_warmup_s = self.time_record_after_warmup_min * 60
        self.iterer_record_after_warmup= np.arange(1, self.time_record_after_warmup_s, 1)
        self.repetitions_record_after_warmup = len(self.iterer_record_after_warmup)


    def initialize_experiment(self) -> TNone:
        f70_status = self.f70_compressor.get_status()
        print(f70_status)

    @rpc
    def cleanup_devices(self) -> TNone:
        """
        Set all devices to normal states
        """
        print("\tCLEANUP - BEGIN")
        # ensure f70 is turned back on
        self.f70_compressor.turn_on()

    """
    MAIN LOGIC
    """
    @rpc
    def run_main(self) -> TNone:
        """
        Run till for duration of experiment or till termination requested
        """

        print('\tTurning OFF F70 Compressor')
        self.f70_compressor.turn_off()
        time_start_quick_warmup = time()
        _loop_iter = 0
        for time_segment in self.iterer_quick_warmup:
            self.check_termination()
            self.record_labrad_values(_loop_iter)
            _loop_iter += 1
            self.update_completion_monitor(
                time_segment,
                self.repetitions_quick_warmup
            )

        time_end_quick_warmup = time()
        print(f'Time of Quick Warmup was : {time_end_quick_warmup - time_start_quick_warmup: .2f} s')
        print('\tTurning ON F70 Compressor')
        self.f70_compressor.turn_on()

        print('\tRecording System Response After Warmup')
        for time_segment in self.iterer_record_after_warmup:
            self.check_termination()
            self.record_labrad_values(_loop_iter)
            _loop_iter += 1
            self.update_completion_monitor(
                time_segment,
                self.repetitions_record_after_warmup
            )

    @rpc
    def record_labrad_values(self,
                             time_s):
        temp_diodes_kelvin = self.lakeshore.get_temperature()

        time_start_recording = time()
        temp_diode_A_kelvin = temp_diodes_kelvin[0]
        temp_diode_B_kelvin = temp_diodes_kelvin[1]
        temp_diode_C_kelvin = temp_diodes_kelvin[2]
        temp_diode_D_kelvin = temp_diodes_kelvin[3]
        pressure_ion_pump_mbar = self.ion_pump.get_pressure()
        time_end_recording = time()
        time_latency_s = time_end_recording - time_start_recording
        sleep(1. - time_latency_s)

        self.update_results(
            time_s,
            temp_diode_A_kelvin,
            temp_diode_B_kelvin,
            temp_diode_C_kelvin,
            temp_diode_D_kelvin,
            pressure_ion_pump_mbar,
        )

    @rpc
    def check_termination(self) -> TNone:
        """
        OVERRIDE base check_termination to ensure devices always clean up
        in case of termination.
        """
        if self.scheduler.check_termination():
            self.cleanup_devices()
            raise TerminationRequested

    @rpc
    def analyze_experiment(self):

        results = self.results
        time_s = self.results[:, 0]
        temp_diode_A_kelvin = results[:,1]
        temp_diode_B_kelvin = results[:,2]
        temp_diode_C_kelvin = results[:,3]
        temp_diode_D_kelvin = results[:, 4]
        pressure_ion_pump_mbar = results[:,5]
        mbar_to_torr = 0.7500062
        pressure_ion_pump_mtorr = (pressure_ion_pump_mbar * mbar_to_torr) / 1e3


        self.set_dataset('pressure_max_ion_pump_mtorr', np.max(pressure_ion_pump_mtorr))

        # format dictionary for applet plotting
        ion_pump_plotting_results = {'x': time_s,
                            'y': pressure_ion_pump_mtorr,
                            'subplot_titles': f'Ion Pump Pressure',
                            'subplot_x_labels': 'Time (s)',
                            'subplot_y_labels': 'Ion Pump Pressure (mTorr)',
                            'rid': self.scheduler.rid,
                             'textbox_strs': [None],
                            }

        self.create_matplotlib_applet(ion_pump_plotting_results,
                                      name=f'Ion Pump Quick Warmup',
                                      group = ['plotting', 'utilities', 'cryovac'],)

        # format dictionary for laekshore applet plotting
        temp_arrs_kevin = np.array([
                                   temp_diode_A_kelvin,
                                   temp_diode_B_kelvin,
                                   temp_diode_C_kelvin,
                                   temp_diode_D_kelvin
                        ])

        time_arrs_s = np.array([
                                time_s,
                                time_s,
                                time_s,
                                time_s,
                                ])

        lakeshore_plotting_results = {'x': time_arrs_s,
                            'y': temp_arrs_kevin,
                            'subplot_titles': f'Lakeshore Temperatures',
                            'subplot_x_labels': 'Time (s)',
                            'subplot_y_labels': 'Diode Temperature (K)',
                            'rid': self.scheduler.rid,
                             'textbox_strs': [None],
                            'legend_labels': ['Diode A', 'Diode B', 'Diode C', 'Diode D'],
                            }

        self.create_matplotlib_applet(lakeshore_plotting_results,
                                      name=f'Lakeshore Temperature Quick Warmup',
                                      group = ['plotting', 'utilities', 'cryovac'],)

        return results

    @rpc(flags={"async"})
    def update_completion_monitor(self,
                                  rep,
                                  segment_length) -> TNone:
        """
        Update the completion monitor to inform user of elapsed time.
        :param rep: recording points
        :param segment_length: length of segment
        """
        self.set_dataset('management.dynamic.completion_pct',
                         round(rep / segment_length * 100., 3),
                         broadcast=True, persist=True, archive=False)

