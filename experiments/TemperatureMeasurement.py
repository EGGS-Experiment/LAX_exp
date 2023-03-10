import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import AbsorptionProbe


class TemperatureMeasurement(LAXExperiment, Experiment):
    """
    Experiment: Temperature Measurement

    Measure the ion temperature by doing a weak probe linescan and fitting the resulting lineshape.
    """

    name = 'Temperature Measurement'

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                            NumberValue(default=100, ndecimals=0, step=1, min=1, max=10000))

        # probe frequency scan
        self.setattr_argument("freq_probe_scan_mhz",                    Scannable(
                                                                            default=RangeScan(85, 135, 10, randomize=True),
                                                                            global_min=85, global_max=135, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=6
                                                                        ))
        # relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('pmt')

        # subsequences
        self.probe_subsequence =                                        AbsorptionProbe(self)

    def prepare_experiment(self):
        # convert probe frequency scans
        self.freq_probe_scan_mhz =                                      np.array([freq_mhz for freq_mhz in self.freq_probe_scan_mhz])
        self.freq_probe_scan_ftw =                                      np.array([hz_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_probe_scan_mhz])

        # get amplitude calibration curve from dataset maanger and interpolate the points
        # interpolation is necessary to allow continuous range of frequency values
        from scipy.interpolate import Akima1DInterpolator
        ampl_calib_points =                                             self.get_dataset('calibration.temperature.asf_calibration_curve_mhz_pct')
        ampl_calib_curve =                                              Akima1DInterpolator(ampl_calib_points[:, 0], ampl_calib_points[:, 1])

        # get calibrated amplitude values
        ampl_probe_scan_pct =                                           np.array(ampl_calib_curve(self.freq_probe_scan_mhz))
        self.ampl_probe_scan_asf =                                      np.array([pct_to_asf(ampl_pct) for ampl_pct in ampl_probe_scan_pct])

        # set up probe waveform config
        self.waveform_probe_scan =                                      np.stack(np.array([self.freq_probe_scan_ftw, self.ampl_probe_scan_asf])).transpose()

    @property
    def results_shape(self):
        return (self.repetitions * len(self.freq_probe_scan_ftw) * 2,
                3)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # record subsequences onto DMA
        self.probe_subsequence.record_dma()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        for trial_num in range(self.repetitions):

            self.core.break_realtime()

            # sweep frequency
            for waveform_params in self.waveform_probe_scan:

                # get waveform parameters
                freq_ftw = waveform_params[0]
                ampl_asf = waveform_params[1]

                # set probe beam frequency
                self.pump.set_mu(freq_ftw, asf=ampl_asf, profile=1)
                self.core.break_realtime()

                # probe absorption at detuning (repump on)
                self.probe_subsequence.run_dma()

                # get counts and update datasets
                counts_actual = self.pmt.fetch_count()
                counts_control = self.pmt.fetch_count()

                # update datasets
                with parallel:
                    self.core.break_realtime()
                    self.update_results(freq_ftw, 1, counts_actual)
                    self.update_results(freq_ftw, 0, counts_control)
