import numpy as np
from artiq.experiment import *

import os
import time
import skimage
from scipy import ndimage
from datetime import datetime
from matplotlib import pyplot as plt

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import Readout
import cv2 as cv


# todo: finish kernel_invariants
# todo: speed things up here and there


class IonLoadAndAramp(LAXExperiment, Experiment):
    """
    Utility: Ion Load and Aramp

    Automatically load ions via resistive oven and eject excess via A-ramp.
    """
    name = 'Ion Load and Aramp'
    BASE_PATH = r"\\eric.physics.ucla.edu\groups\motion\Data"

    kernel_invariants = {
        "pmt_sample_num", "pmt_dark_threshold_counts",
        "att_397_mu", "att_866_mu", "att_854_mu",
        "time_runtime_max_s", "start_time_s",
        "IMAGE_HEIGHT", "IMAGE_WIDTH", "image_region", "data_path",
        "time_aramp_pulse_s"
    }

    def build_experiment(self):
        # general arguments
        self.setattr_argument('desired_num_of_ions', NumberValue(default=2, min=1, max=10, precision=0, step=1))

        # starting trap arguments
        self.setattr_argument('start_east_endcap_voltage',  NumberValue(default=18.8, precision=1, step=0.1, min=0., max=300.),
                                                            group='Starting Trap Parameters')
        self.setattr_argument('start_west_endcap_voltage',  NumberValue(default=24., precision=1, step=0.1, min=0., max=300.),
                                                            group='Starting Trap Parameters')

        # ending trap arguments
        self.setattr_argument('end_east_endcap_voltage',     NumberValue(default=205., precision=1, step=0.1, min=0., max=400.),
                                                                group='Ending Trap Parameters')
        self.setattr_argument('end_west_endcap_voltage',     NumberValue(default=308., precision=1, step=0.1, min=0., max=400.),
                                                                group='Ending Trap Parameters')
        self.setattr_argument('end_v_shim_voltage',          NumberValue(default=67.6, precision=1, step=0.1, min=0., max=150.),
                                                                group='Ending Trap Parameters')
        self.setattr_argument('end_h_shim_voltage',          NumberValue(default=48.3, precision=1, step=0.1, min=0., max=150.),
                                                                group='Ending Trap Parameters')
        self.setattr_argument('end_aramp_voltage',          NumberValue(default=2.8, precision=1, step=0.1, min=0., max=50.),
                                                                group='Ending Trap Parameters')

        # aramping parameters
        self.setattr_argument("enable_aramp",               BooleanValue(default=False), group='A-Ramp Ejection')
        self.setattr_argument("aramp_ions_voltage_list",    Scannable(
                                                                    default=[
                                                                        RangeScan(18, 24, 20, randomize=True),
                                                                        ExplicitScan([19, 20, 21, 22, 23, 24]),
                                                                    ],
                                                                    global_min=0.0, global_max=30.0, global_step=1,
                                                                    unit="V", scale=1, precision=2
                                                                ), group='A-Ramp Ejection')

        # oven configuration
        self.setattr_argument("oven_voltage",   NumberValue(default=1.25, max = 2., precision=2, step=0.01, unit="V"),
        self.setattr_argument("oven_voltage",   NumberValue(default=1.25, max = 2., precision=2, step=0.01, unit="V"),
                                                        group='Oven Settings')
        self.setattr_argument("oven_current",   NumberValue(default=3.25, max=4., precision=2, step=0.01, unit="A"),
                                                        group='Oven Settings')

        # image region parameters: MAX (450,450) TO PREVENT LASER SCATTER OFF ELECTRODES FROM CONFUSING ANALYSIS
        self.setattr_argument('image_width_pixels',     NumberValue(default=400, min=100, max=450, step=50, scale=1, precision=0), group='Camera')
        self.setattr_argument('image_height_pixels',    NumberValue(default=400, min=100, max=450, step=50, scale=1, precision=0), group='Camera')
        self.setattr_argument('horizontal_binning',     NumberValue(default=1, min=1, max=5, step=1, scale=1, precision=0), group='Camera')
        self.setattr_argument('vertical_binning',       NumberValue(default=1, min=1, max=5, step=1, scale=1, precision=0), group='Camera')

        # relevant devices - sinara
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('scheduler')
        self.readout_subsequence = Readout(self)

        # relevant devices - labrad
        self.setattr_device('shutters')
        self.setattr_device('oven')
        self.setattr_device('aperture')
        self.setattr_device('trap_dc')
        self.setattr_device('camera')
        self.setattr_device('flipper')

    def prepare_experiment(self):
        """
        Prepare experimental values and precompute/preallocate
        to reduce kernel overheads.
        """
        '''HARDWARE VALUES'''
        self.pmt_sample_num = 30
        self.pmt_dark_threshold_counts = 15
        self.time_runtime_max_s = 900.

        self.start_time_s = np.int64(0)

        '''CAMERA SETUP'''
        self.IMAGE_HEIGHT = 512
        self.IMAGE_WIDTH =  512
        border_x, border_y = ((self.IMAGE_WIDTH - self.image_width_pixels) / 2,
                              (self.IMAGE_HEIGHT - self.image_height_pixels) / 2)
        start_x, start_y =  (int(border_x), int(border_y))
        end_x, end_y =      (int(self.IMAGE_WIDTH - border_x) - 1,
                             int(self.IMAGE_WIDTH - border_y) - 1)
        self.image_region = (self.horizontal_binning, self.vertical_binning,
                             start_x, end_x, start_y, end_y)

        '''BEAM SETUP'''
        # set attenuations (per Josh's request)
        self.att_397_mu = att_to_mu(14 * dB)
        self.att_854_mu = att_to_mu(14 * dB)
        self.att_866_mu = att_to_mu(14 * dB)

        '''A-RAMP SETUP'''
        self.time_aramp_pulse_s = 2.

        '''FILEPATHS'''
        # todo: migrate data storage to a dataset so we can re-process images later
        year_month =        datetime.today().strftime('%Y-%m')
        year_month_day =    datetime.today().strftime('%Y-%m-%d')
        self.data_path = os.path.join(self.BASE_PATH, year_month, year_month_day)
        # ensure data folder exists; create it if it doesn't exist
        os.makedirs(self.data_path, exist_ok=True)

    @property
    def results_shape(self):
        return (2, 2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        self.core.break_realtime()

        # store attenuations to prevent overriding
        self.pump.beam.cpld.get_att_mu()
        self.core.break_realtime()

        # set readout profile for beams
        self.pump.readout()
        self.core.break_realtime()
        self.pump.set_att_mu(self.att_397_mu)
        self.repump_qubit.set_att_mu(self.att_854_mu)
        self.repump_cooling.set_att_mu(self.att_866_mu)
        self.core.break_realtime()
        # turn on lasers
        self.pump.on()
        self.repump_qubit.on()
        self.repump_cooling.on()
        self.core.break_realtime()

        # deterministically set flipper to camera
        self.set_flipper_to_camera()

        # synchronize timeline
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

    @rpc
    def initialize_labrad_devices(self) -> TNone:
        """
        Initialize labrad devices via RPC.
        """
        try:
            '''SET UP CAMERA'''
            # set camera region of interest and exposure time
            self.camera.set_image_region(self.image_region)

            '''SET UP TRAP'''
            # set endcaps to loading voltages
            self.trap_dc.set_east_endcap_voltage(self.start_east_endcap_voltage)
            self.trap_dc.set_west_endcap_voltage(self.start_west_endcap_voltage)
            # turn on endcap channels and ensure others are off
            self.trap_dc.east_endcap_toggle(True)
            self.trap_dc.west_endcap_toggle(True)
            self.trap_dc.h_shim_toggle(False)
            self.trap_dc.v_shim_toggle(False)
            self.trap_dc.aramp_toggle(False)

            '''SET UP LOADING LASERS'''
            # open 397nm aperture
            self.aperture.open_aperture()
            # open shutters
            self.shutters.toggle_377_shutter(True)
            self.shutters.toggle_423_shutter(True)

            '''START OVEN'''
            # turn on the oven
            self.oven.set_oven_voltage(self.oven_voltage)
            self.oven.set_oven_current(self.oven_current)
            self.oven.toggle(True)

        except Exception as e:
            print(repr(e))
            self.cleanup_devices()

    @rpc
    def run_main(self) -> TNone:
        """
        Run till ion is loaded or timeout.
        """
        # get start time to check if we exceed max time
        self.start_time_s = time.time()

        # run loading loop until we load desired_num_of_ions
        num_ions = 0

        try:
            while (num_ions != self.desired_num_of_ions) and (num_ions != -1):
                self.check_termination()

                # load ions if below desired count
                if num_ions < self.desired_num_of_ions:
                    self.initialize_labrad_devices()
                    num_ions = self.load_ion()

                # eject excess ions via A-ramping (if enabled)
                elif self.enable_aramp and (num_ions > self.desired_num_of_ions):
                    self.cleanup_devices()
                    num_ions = self.aramp_ions()

                # otherwise, simply stop execution
                elif (not self.enable_aramp) and (num_ions > self.desired_num_of_ions):
                    print("\tTOO MANY IONS LOADED - STOPPING HERE.")
                    break

            print("\tLOADING COMPLETE - CLEANING UP.")

        except Exception as e:
            print("Error - stopping execution; cleaning up.")
            print(repr(e))

        '''CLEAN UP'''
        self.cleanup_devices()

    @rpc
    def load_ion(self) -> TInt32:
        """
        Function for loading an ion.
        Returns:
            TInt32: number of ions in trap
        """
        print("\tLOAD ION - BEGIN:")

        # define some variables for later use
        idx = 0
        ion_spottings = 0

        # run loop while we don't see ions
        while True:
            # extract number of ions on camera
            num_ions = self.process_image('pre_aramp_original.png', 'pre_aramp_maipulated.png')

            # periodically check if we've reached a stop condition
            # i.e. max_time reached or termination_requested
            if idx % 5 == 0:
                if (time.time() - self.start_time_s) > self.time_runtime_max_s:
                    print("\t\tPROBLEM: TOOK OVER 15 MIN TO LOAD - ENDING PROGRAM")
                    return -1
                else:
                    self.check_termination()

            # check to see if we have loaded enough ions
            if num_ions >= self.desired_num_of_ions:
                ion_spottings += 1
                # ensure camera sees ion in 3 consecutive images to prevent singular false positive
                if ion_spottings >= 3:
                    print("\t\t{:d} ION(s) LOADED".format(num_ions))
                    return num_ions
            else:
                ion_spottings = 0  # reset if image analysis shows no ions in trap
            idx += 1

        return 0

    @rpc
    def aramp_ions(self) -> TInt32:
        """
        Pulse A-ramp to get rid of excess ions.
        Returns:
            TInt32: Number of ions.
        """
        # todo: allow for some error conditions
        print("\tARAMP EJECTION - START:")
        for aramp_voltage in self.aramp_ions_voltage_list:
            self.check_termination()

            # pulse A-ramp for given period
            print("\t\tARAMPING @ VOLTAGE: {:.1f} V".format(aramp_voltage))
            self.trap_dc.set_aramp_voltage(aramp_voltage)
            time.sleep(self.time_aramp_pulse_s)
            self.trap_dc.set_aramp_voltage(self.end_aramp_voltage)
            time.sleep(self.time_aramp_pulse_s)

            # get number of ions
            num_ions = self.process_image('post_aramp_original.png', 'post_aramp_manipulated.png')
            print("\t\tIONS REMAINING: {:d}".format(num_ions))
            if num_ions <= self.desired_num_of_ions:
                return num_ions

        return 0

    @rpc
    def cleanup_devices(self) -> TNone:
        """
        Set all devices to states as if ion was loaded.
        All but trap electrodes set to original state --- trap electrodes set to final trapping potential.
        """
        # turn off oven
        self.oven.set_oven_voltage(0.)
        self.oven.set_oven_current(0.)
        self.oven.toggle(False)

        # close shutters
        self.shutters.toggle_377_shutter(False)
        self.shutters.toggle_423_shutter(False)

        # set trap parameters as if ion was loaded
        self.trap_dc.set_east_endcap_voltage(self.start_east_endcap_voltage)
        self.trap_dc.set_west_endcap_voltage(self.start_west_endcap_voltage)
        self.trap_dc.set_h_shim_voltage(self.end_h_shim_voltage)
        self.trap_dc.set_v_shim_voltage(self.end_v_shim_voltage)
        self.trap_dc.set_aramp_voltage(self.end_aramp_voltage)

        # turn on the endcap channels
        self.trap_dc.h_shim_toggle(True)
        self.trap_dc.v_shim_toggle(True)
        self.trap_dc.aramp_toggle(True)

        # ramp endcaps to end values
        self.trap_dc.ramp_both_endcaps([self.end_east_endcap_voltage, self.end_west_endcap_voltage],
                                       [100., 100.])
        # note: we add sleep and set voltages AGAIN to reflect updates on GUIs
        time.sleep(2)
        self.trap_dc.set_east_endcap_voltage(self.end_east_endcap_voltage)
        self.trap_dc.set_west_endcap_voltage(self.end_west_endcap_voltage)

        # close aperture
        # note: do this last to keep ion happy during voltage adjustments
        self.aperture.close_aperture()

    @rpc
    def process_image(self, filepath1: TStr=None, filepath2: TStr=None) -> TInt32:
        """
        Process image data from camera and extract the number of ions in the trap.
        Arguments:
            filepath1   (string)    : save path for the original (no post-processing) image.
            filepath2   (string)    : save path for the processed image.
        Returns:
            TInt32: Number of ions in the trap
        """
        if filepath1 is None:   filepath1 = "original.png"
        if filepath2 is None:   filepath2 = "manipulated.png"

        # get camera data and format into image
        image_arr = self.camera.get_most_recent_image()
        data = np.reshape(image_arr, (self.image_width_pixels, self.image_height_pixels))



        # create camera image (raw)
        plt.figure(1)
        plt.title("Original Image")
        plt.imshow(data)
        plt.savefig(os.path.join(self.data_path, filepath1))

        # # conduct binary thresholding on camera image
        # upper_percentile = np.percentile(data, 99.97)
        # data[data < upper_percentile] =     0
        # data[data >= upper_percentile] =    1
        #
        # # extract number of ions in image
        # kernel1 = np.ones((2, 2))
        # data = ndimage.binary_erosion(data, kernel1, iterations=2)
        # data[data > 0] = 1
        # labels = skimage.measure.label(data)

        data = ((data - np.min(data)) / (np.max(data) - np.min(data))) * 255
        data = np.uint8(data)
        data = data * (data > np.quantile(data, 0.99))
        data = cv.medianBlur(data, 5)

        circles = cv.HoughCircles(data, cv.HOUGH_GRADIENT, 1, 5,
                                  param1=10, param2=5, minRadius=2, maxRadius=5)

        if circles is not None:
            num_circles = len(circles[0])
            circles = np.uint16(np.around(circles))
            for i in circles[0, :]:
                # draw the outer circle
                cv.circle(data, (i[0], i[1]), i[2], (255, 255, 255), 1)
                # draw the center of the circle
                cv.circle(data, (i[0], i[1]), 2, (255, 255, 255), 1)
        else:
            num_ions = 0

        # create camera image (processed)
        plt.figure(2)
        plt.imshow(data)
        plt.title("Manipulated Image")
        plt.savefig(os.path.join(self.data_path, filepath2))
        return num_ions

    @kernel(flags={"fast-math"})
    def set_flipper_to_camera(self) -> TNone:
        """
        Set the flipper to send light to the camera.
        """
        # store PMT count events
        for i in range(self.pmt_sample_num):
            self.readout_subsequence.run()
            delay_mu(100)
        self.core.break_realtime()

        # retrieve counts from PMT
        counts = 0
        for i in range(self.pmt_sample_num):
            counts += self.readout_subsequence.fetch_count()
            delay_mu(100)

        # if PMT counts are below the dark count threshold flip again
        if counts > (self.pmt_dark_threshold_counts * self.pmt_sample_num):
            self.flipper.flip()
        self.core.break_realtime()

    @rpc
    def check_termination(self) -> TNone:
        """
        OVERRIDE base check_termination to ensure devices are always cleaned up
        in case of termination.
        """
        if self.scheduler.check_termination():
            self.cleanup_devices()
            raise TerminationRequested
