import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
import LAX_exp.experiments.diagnostics.LaserScan as LaserScan
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS
from skimage.transform import hough_circle, hough_circle_peaks
from time import time, sleep
from matplotlib.pyplot import imsave
from datetime import datetime
import os

class AxialFreqHysteresis(LaserScan.LaserScan):
    """
    Experiment: Axial Freq Hysteresis

    Does a 729nm laser scan; resets the ion(s) every shot.
    Supports sine-squared pulse shaping.
    Scans DC output for axial trapping potential to determine hystersis
    """
    name = 'Axial Freq Hysteresis'
    # note: no need to specify kernel invariants since parent specifies all of them for us

    def build_experiment(self):
        # call parent build
        super().build_experiment()

        self.setattr_argument('dc_voltage_endcap_scan', Scannable([
            ExplicitScan([-1,0,1]),
            RangeScan(-1,1, 10, randomize=False),
            CenterScan(0,5,1, randomize=False)
        ],
            global_max=400,
            global_min=0, global_step=0.01,
            unit = 'V', precision=2
        ), group='voltage_scan',
        tooltip='relative change in voltage to BOTH endcaps'
                'NOTE BENE: these values are how much the endcap voltages are CHANGED by,'
                'NOT THE ABSOLUTE voltages')

        self.setattr_argument("num_voltage_sweeps", NumberValue(1,
                                                                step=1, precision=0,
                                                                max = 100, min=1),
                              group='voltage_scan',
                              tooltip='Number of voltage sweeps to perform when determining if there is hystersis')

        # camera features for calibrating scaling of DC voltages
        self.setattr_argument('image_center_x',
                              NumberValue(default=256, min=1, max=512, step=50, precision=0, scale=1, unit="pixels"),
                              group='Camera',
                              tooltip="Image center position (x).")
        self.setattr_argument('image_center_y',
                              NumberValue(default=256, min=1, max=512, step=50, precision=0, scale=1, unit="pixels"),
                              group='Camera',
                              tooltip="Image center position (y).")
        self.setattr_argument('image_width_pixels',
                              NumberValue(default=350, min=100, max=450, step=50, precision=0, scale=1, unit="pixels"),
                              group='Camera',
                              tooltip="Square width of the total 512x512 image region on camera.\n"
                                      "Defined relative to (image_center_x, image_center_y).\n"
                                      "Note: too large an image width (~450 pixels) may contain scatter and confuse the analysis.")
        self.setattr_argument('horizontal_binning',
                              NumberValue(default=1, min=1, max=5, step=1, precision=0, scale=1, unit="bins"),
                              group='Camera',
                              tooltip="Horizontal pixel bin size to set on camera.\n"
                                      "Note: this is software binning, as opposed to hardware binning.")
        self.setattr_argument('vertical_binning',
                              NumberValue(default=1, min=1, max=5, step=1, precision=0, scale=1, unit="bins"),
                              group='Camera',
                              tooltip="Vertical pixel bin size to set on camera.\n"
                                      "Note: this is software binning, as opposed to hardware binning.")

        # laser scan multi - arguments
        self.setattr_device('trap_dc')
        self.setattr_device('camera')
        self.setattr_device('flipper')

    def prepare_experiment(self):
        super().prepare_experiment()

        self.west_endcap_voltage = self.trap_dc.get_west_endcap_voltage()
        self.east_endcap_voltage = self.trap_dc.get_east_endcap_voltage()
        self.dc_voltage_endcap_scan = list(self.dc_voltage_endcap_scan)

        '''
        CAMERA SETUP
        '''
        # get image height and width from camera arguments (info_detector_dimensions)
        IMAGE_WIDTH, IMAGE_HEIGHT = self.camera.acquire_detector_dimensions()

        # process specified image region in terms of width and position
        start_x, start_y = (
            max(1, self.image_center_x - self.image_width_pixels / 2),
            max(1, self.image_center_y - self.image_width_pixels / 2)
        )
        end_x, end_y = (
            min(IMAGE_WIDTH, self.image_center_x + self.image_width_pixels/2 - 1),
            min(IMAGE_HEIGHT, self.image_center_y + self.image_width_pixels/2 - 1)
        )
        self.image_region = (self.horizontal_binning, self.vertical_binning,
                             int(start_x), int(end_x), int(start_y), int(end_y))
        self.set_dataset("camera_image_region", self.image_region, broadcast=False)

        try:
            '''SET UP CAMERA'''
            # set camera region of interest and exposure time
            self.camera.set_image_region(self.image_region)

        except Exception as e:
            print("Error during initialize_labrad_devices: {:}".format(repr(e)))

        self.path_image_save = r"\\eric.physics.ucla.edu\groups\motion\Data"
        year_month =        datetime.today().strftime('%Y-%m')
        year_month_day =    datetime.today().strftime('%Y-%m-%d')
        self.data_path = os.path.join(self.path_image_save, year_month, year_month_day)
        os.makedirs(self.data_path, exist_ok=True)

    @property
    def results_shape(self):
        return (self.repetitions * self.num_voltage_sweeps * len(self.dc_voltage_endcap_scan) * len(self.config_experiment_list),
                4)

    @kernel(flags={'fast-math'})
    def run_main(self) -> TNone:

        self.core.break_realtime()
        self.core.wait_until_mu(now_mu())
        self.locate_ion()
        self.core.break_realtime()

        for num in range(self.num_voltage_sweeps):
            for voltage in self.dc_voltage_endcap_scan:

                # consume all slack to ensure all events have been submitted
                self.core.break_realtime()
                self.core.wait_until_mu(now_mu())

                self.trap_dc.set_east_endcap_voltage(self.east_endcap_voltage + voltage)
                self.trap_dc.set_west_endcap_voltage(self.west_endcap_voltage + voltage)

                # add slack back after rpc call
                self.core.break_realtime()

                # let voltages settle
                delay_mu(np.int64(2e9))

                for trial_num in range(self.repetitions):
                    for config_vals in self.config_experiment_list:

                        ### PREPARE & CONFIGURE ###
                        # extract values from config list
                        freq_ftw = int32(config_vals[0])
                        time_holdoff_mu = config_vals[1]

                        # prepare relevant beams
                        self.core.break_realtime()
                        if self.enable_pulseshaping:
                            self.qubit.set_ftw(freq_ftw)
                        else:
                            self.qubit.set_mu(freq_ftw, asf=self.ampl_qubit_asf,
                                              profile=self.profile_729_readout,
                                              phase_mode=PHASE_MODE_CONTINUOUS)
                        delay_mu(10000)

                        # wait for linetrigger
                        if self.enable_linetrigger:
                            self.trigger_line.trigger(self.trigger_line.time_timeout_mu, time_holdoff_mu)

                        ### MAIN SHOT ###
                        # initialize ion in S-1/2 state
                        self.initialize_subsequence.run_dma()

                        # fire spectroscopy pulse
                        if self.enable_pulseshaping:
                            self.pulseshape_subsequence.run()
                        else:
                            self.qubit.on()
                            delay_mu(self.time_qubit_mu)
                            self.qubit.off()

                        # read out counts & clean up loop
                        self.readout_subsequence.run_dma()
                        self.rescue_subsequence.resuscitate()
                        self.initialize_subsequence.slack_rescue()
                        counts = self.readout_subsequence.fetch_count()

                        # store results in dataset
                        self.rescue_subsequence.detect_death(counts)
                        self.update_results(freq_ftw, counts, time_holdoff_mu, voltage)

                        # rescue ion as needed & support graceful termination
                        self.core.break_realtime()
                        self.rescue_subsequence.run(trial_num)
                        self.check_termination()

    def analyze_experiment(self):
        pass

    def cleanup_experiment(self) -> TNone:

        self.trap_dc.set_east_endcap_voltage(self.east_endcap_voltage)
        self.trap_dc.set_west_endcap_voltage(self.west_endcap_voltage)


    @rpc
    def locate_ion(self) -> TNone:
        """
        Process image data from camera and extract location of ion
        Returns:
            TInt32: Number of ions in the trap
        """

        # get camera data and reshape into image
        # self.camera.acquire_single_image()
        image_arr = self.camera.get_most_recent_image()
        data = np.reshape(image_arr, (self.image_width_pixels, self.image_width_pixels))
        imsave(os.path.join(self.data_path, "hystersis.jpg"), data)

        # threshold & rescale data
        # todo: set 1000 as some parameter for min scatter value
        data *= data > 1000
        data = np.uint8(((data - np.min(data)) / (np.max(data) - np.min(data))) * 255)
        # use only upper 1% quantile of data
        data *= data > np.quantile(data, 0.99)
        imsave(os.path.join(self.data_path, "hystersis2.jpg"), data)

        # extract ion positions
        guess_radii = np.arange(1, 8)
        circles = hough_circle(data, guess_radii)
        accums, cxs, cys, radii = hough_circle_peaks(
            circles, guess_radii, min_xdistance=1, min_ydistance=1, threshold=0.95)

        # create unique list of cx, cy coordinates
        unique_locs = set(tuple(
            (cx, cys[idx])
            for idx, cx in enumerate(cxs)
        ))

        print(unique_locs)