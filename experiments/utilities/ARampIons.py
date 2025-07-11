import numpy as np
from artiq.experiment import *

import os
import time
from datetime import datetime
from matplotlib import pyplot as plt

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import Readout

import skimage
from skimage import color
from skimage.transform import hough_circle, hough_circle_peaks


class ARampEjection(LAXExperiment, Experiment):
    """
    Utility: A-Ramp Ejection

    Eject excess ions in the trap by pulsing the voltage on the A-ramp.
    """
    name = 'Aramp'
    BASE_PATH = r"\\eric.physics.ucla.edu\groups\motion\Data"

    kernel_invariants = {
        "pmt_sample_num", "pmt_dark_threshold_counts",
        "att_397_mu", "att_866_mu", "att_854_mu",
        "IMAGE_HEIGHT", "IMAGE_WIDTH", "image_region", "data_path", "time_aramp_pulse_s"
    }

    def build_experiment(self):
        # general arguments
        self.setattr_argument('desired_num_of_ions', NumberValue(default=1, min=1, max=10, precision=0, step=1))


        # ending trap arguments
        self.setattr_argument('east_endcap_voltage',     NumberValue(default=205., precision=1, step=0.1, min=0., max=400.),
                                                                group='Ending Trap Parameters')
        self.setattr_argument('west_endcap_voltage',     NumberValue(default=308., precision=1, step=0.1, min=0., max=400.),
                                                                group='Ending Trap Parameters')
        self.setattr_argument('v_shim_voltage',          NumberValue(default=70.8, precision=1, step=0.1, min=0., max=150.),
                                                                group='Ending Trap Parameters')
        self.setattr_argument('h_shim_voltage',          NumberValue(default=50.5, precision=1, step=0.1, min=0., max=150.),
                                                                group='Ending Trap Parameters')
        self.setattr_argument('final_aramp_voltage',          NumberValue(default=2.8, precision=1, step=0.1, min=0., max=50.),
                                                                group='Ending Trap Parameters')

        # image region parameters: MAX (450,450) TO PREVENT LASER SCATTER OFF ELECTRODES FROM CONFUSING ANALYSIS
        self.setattr_argument('image_width_pixels',     NumberValue(default=400, min=100, max=450, step=50, scale=1, precision=0), group='Camera')
        self.setattr_argument('image_height_pixels',    NumberValue(default=400, min=100, max=450, step=50, scale=1, precision=0), group='Camera')
        self.setattr_argument('horizontal_binning',     NumberValue(default=1, min=1, max=5, step=1, scale=1, precision=0), group='Camera')
        self.setattr_argument('vertical_binning',       NumberValue(default=1, min=1, max=5, step=1, scale=1, precision=0), group='Camera')

        # aramping parameters
        self.setattr_argument("aramp_ions_voltage_list",    Scannable(
                                                                    default=[
                                                                        ExplicitScan([14, 14.5, 15, 15.5, 16]),
                                                                        RangeScan(18, 24, 20, randomize=True),
                                                                    ],
                                                                    global_min=0.0, global_max=30.0, global_step=1,
                                                                    unit="V", scale=1, precision=2
                                                                ), group='A-Ramp Ejection')

        # relevant devices - sinara
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('scheduler')
        self.readout_subsequence = Readout(self)

        # relevant devices - labrad
        self.setattr_device('trap_dc')
        self.setattr_device('camera')
        self.setattr_device('flipper')
        self.setattr_device('aperture')

    def prepare_experiment(self):
        """
        Prepare experimental values and precompute/preallocate
        to reduce kernel overheads.
        """
        '''HARDWARE VALUES'''
        self.pmt_sample_num = 30
        self.pmt_dark_threshold_counts = 15

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

    @property
    def results_shape(self):
        return (2, 2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        """
        Hardware initialization immediately before run().
        """
        # store attenuations to prevent overriding
        self.pump.beam.cpld.get_att_mu()
        self.core.break_realtime()

        # set readout profile for beams
        self.pump.readout()
        self.pump.set_att_mu(self.att_397_mu)
        delay_mu(8000)
        self.repump_qubit.set_att_mu(self.att_854_mu)
        delay_mu(500)
        self.repump_cooling.set_att_mu(self.att_866_mu)
        delay_mu(8000)
        # turn on lasers
        self.pump.on()
        self.repump_qubit.on()
        self.repump_cooling.on()
        delay_mu(8000)

        # deterministically set flipper to camera
        self.set_flipper_to_camera()

        '''SET UP CAMERA'''
        # set camera region of interest and exposure time
        self.camera.set_image_region(self.image_region)
        self.core.break_realtime()

        '''SET UP TRAP'''
        # ramp endcaps in case not already ramped
        self.trap_dc.ramp_both_endcaps([self.east_endcap_voltage, self.west_endcap_voltage],
                                       [100., 100.])
        # note: we add sleep and set voltages AGAIN to reflect updates on GUIs
        time.sleep(3)
        self.trap_dc.set_east_endcap_voltage(self.east_endcap_voltage)
        self.trap_dc.set_west_endcap_voltage(self.west_endcap_voltage)
        self.trap_dc.set_h_shim_voltage(self.h_shim_voltage)
        self.trap_dc.set_v_shim_voltage(self.v_shim_voltage)
        # turn on voltages
        self.trap_dc.east_endcap_toggle(True)
        self.trap_dc.west_endcap_toggle(True)
        self.trap_dc.h_shim_toggle(False)
        self.trap_dc.v_shim_toggle(False)
        self.trap_dc.aramp_toggle(True)
        self.core.break_realtime()

        # synchronize timeline
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

    def run_main(self) -> TNone:
        """
        Run till ion is loaded or timeout.
        """
        self.aramp_ions()
        self.cleanup_devices()

    @rpc
    def aramp_ions(self) -> TNone:
        """
        Pulse A-ramp to get rid of excess ions.
        """
        # todo: allow for some error conditions
        print("STARTING ARAMP PROCEDURE")

        self.aperture.open_aperture()
        for aramp_voltage in self.aramp_ions_voltage_list:

            # get number of ions
            num_ions = self.process_image('post_aramp_original.png', 'post_aramp_manipulated.png')
            if num_ions <= self.desired_num_of_ions:
                break

            # pulse A-ramp for given period
            print(f"ARAMPING AT VOLTAGE {np.round(aramp_voltage,2)}")
            self.trap_dc.set_aramp_voltage(aramp_voltage)
            time.sleep(self.time_aramp_pulse_s)
            self.trap_dc.set_aramp_voltage(self.final_aramp_voltage)
            time.sleep(self.time_aramp_pulse_s)

            # check termination
            self.check_termination()

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

        # todo: set 1000 as some parameter for min scattering value
        data = data * (data > 1000)
        if np.max(data) > 0:
            data = ((data - np.min(data)) / (np.max(data) - np.min(data))) * 255
        data = np.uint8(data)
        data = data * (data > np.quantile(data, 0.99))

        guess_radii = np.arange(1, 8)
        circles = hough_circle(data, guess_radii)
        accums, cx, cy, radii = hough_circle_peaks(circles, guess_radii, min_xdistance=1, min_ydistance=1,
                                                   threshold=0.95)
        num_ions = len(cx)

        data = color.gray2rgb(data)
        for center_x, center_y, radius in zip(cx, cy, radii):
            circy, circx = skimage.draw.circle_perimeter(center_x, center_y, radius, shape=data.shape)
            data[circx, circy] = (220, 20, 20)

        # create camera image (processed)
        plt.figure(2)
        plt.imshow(data)
        plt.title("Manipulated Image")
        plt.savefig(os.path.join(self.data_path, filepath2))
        return num_ions

    @rpc
    def cleanup_devices(self) -> TNone:
        """
        Set all devices to states as if ion was loaded.
        All but trap electrodes set to original state --- trap electrodes set to final trapping potential
        """
        # set final trap parameters
        self.trap_dc.set_aramp_voltage(self.final_aramp_voltage)

        # turn on the endcap channels
        self.trap_dc.h_shim_toggle(True)
        self.trap_dc.v_shim_toggle(True)
        self.trap_dc.aramp_toggle(True)

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
