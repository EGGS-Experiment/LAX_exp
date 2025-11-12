import numpy as np
from artiq.experiment import *
from matplotlib.pyplot import imsave

import os
from time import time, sleep
from datetime import datetime

from LAX_exp.language import *
from LAX_exp.system.subsequences import Readout
from skimage.transform import hough_circle, hough_circle_peaks
# todo: save image before anything
# todo: save image after ready to load
# todo: save image immediately before ion detected (store running holder of old data)
# todo: save image after ion detected but not yet cleaned up
# todo: save image after ion detected and after cleanup


class IonLoadAndAramp(LAXExperiment, Experiment):
    """
    Utility: Ion Load and Aramp

    Automatically load ions via resistive oven and eject excess via A-ramp.
    """
    name = 'Ion Load and Aramp'

    kernel_invariants = {
        # hardware values

        # magic numbers
        "pmt_sample_num", "pmt_dark_threshold_counts", "pmt_flip_threshold_std",
        "cam_sample_num", "cam_flip_threshold_avg",
        "time_runtime_max_s", "start_time_s", "time_aramp_pulse_s",

        # camera configuration
        "image_region", "path_image_save", "data_path",
    }

    def build_experiment(self):
        # general arguments
        self.setattr_argument('desired_num_of_ions', NumberValue(default=1, min=1, max=10, precision=0, step=1),
                              tooltip="Number of ions to load.\n"
                                      "Ion number is determined by camera image recognition.")

        # starting trap arguments
        self.setattr_argument('start_east_endcap_voltage',  NumberValue(default=16, precision=1, step=0.1, min=0., max=400., scale=1., unit="V"),
                              group='Start Trap Params',
                              tooltip="East endcap voltage to use for loading.\n"
                                      "Lower voltages reduce the trap depth, making it easier to load.\n"
                                      "However, too low a trap depth is unstable, hindering loading.")
        self.setattr_argument('start_west_endcap_voltage',  NumberValue(default=35., precision=1, step=0.1, min=0., max=400., scale=1., unit="V"),
                              group='Start Trap Params',
                              tooltip="West endcap voltage to use for loading.\n"
                                      "Lower voltages reduce the trap depth, making it easier to load.\n"
                                      "However, too low a trap depth is unstable, hindering loading.")

        # ending trap arguments
        self.setattr_argument('end_east_endcap_voltage',    NumberValue(default=177., precision=1, step=0.1, min=0., max=400., scale=1., unit="V"),
                              group='End Trap Params')
        self.setattr_argument('end_west_endcap_voltage',    NumberValue(default=353., precision=1, step=0.1, min=0., max=400., scale=1., unit="V"),
                              group='End Trap Params')
        self.setattr_argument('end_v_shim_voltage',         NumberValue(default=62.8, precision=1, step=0.1, min=0., max=150., scale=1., unit="V"),
                              group='End Trap Params')
        self.setattr_argument('end_h_shim_voltage',         NumberValue(default=78.3, precision=1, step=0.1, min=0., max=150., scale=1., unit="V"),
                              group='End Trap Params')
        self.setattr_argument('end_aramp_voltage',          NumberValue(default=2.8, precision=1, step=0.1, min=0., max=50., scale=1., unit="V"),
                              group='End Trap Params')

        # aramping parameters
        self.setattr_argument("enable_aramp",               BooleanValue(default=False),
                              group='A-Ramp',
                              tooltip="Enable automatic A-Ramp ejection of excess ions.\n"
                                      "Ion number and ejection are determined by camera image recognition.")
        self.setattr_argument("aramp_ions_voltage_list",    Scannable(
                                                                    default=[
                                                                        RangeScan(16, 17.5, 20, randomize=False),
                                                                        ExplicitScan([19, 20, 21, 22, 23, 24]),
                                                                    ],
                                                                    global_min=0.0, global_max=30.0, global_step=1,
                                                                    unit="V", scale=1, precision=2
                                                                ),
                              group='A-Ramp',
                              tooltip="The list of voltages to scan when ejecting excess ions.")

        # oven configuration
        self.setattr_argument("oven_voltage",   NumberValue(default=1.50, min=0., max=2., precision=2, step=0.1, unit="V", scale=1.),
                              group='Oven',
                              tooltip="todo: document")
        self.setattr_argument("oven_current",   NumberValue(default=3.25, min=0., max=4., precision=2, step=0.1, unit="A", scale=1.),
                              group='Oven',
                              tooltip="todo: document")

        # image region parameters: MAX (450,450) TO PREVENT LASER SCATTER OFF ELECTRODES FROM CONFUSING ANALYSIS
        self.setattr_argument("set_to_pmt_after_loading", BooleanValue(True),
                              group='Camera',
                              tooltip="Set the flipper to the PMT after loading.\n"
                                      "Useful for pipelining experiment operation.")
        self.setattr_argument('image_center_x', NumberValue(default=256, min=1, max=512, step=50, precision=0, scale=1, unit="pixels"),
                              group='Camera',
                              tooltip="Image center position (x).")
        self.setattr_argument('image_center_y', NumberValue(default=256, min=1, max=512, step=50, precision=0, scale=1, unit="pixels"),
                              group='Camera',
                              tooltip="Image center position (y).")
        self.setattr_argument('image_width_pixels', NumberValue(default=350, min=100, max=450, step=50, precision=0, scale=1, unit="pixels"),
                              group='Camera',
                              tooltip="Square width of the total 512x512 image region on camera.\n"
                                      "Defined relative to (image_center_x, image_center_y).\n"
                                      "Note: too large an image width (~450 pixels) may contain scatter and confuse the analysis.")
        self.setattr_argument('horizontal_binning', NumberValue(default=1, min=1, max=5, step=1, precision=0, scale=1, unit="bins"),
                              group='Camera',
                              tooltip="Horizontal pixel bin size to set on camera.\n"
                                      "Note: this is software binning, as opposed to hardware binning.")
        self.setattr_argument('vertical_binning',   NumberValue(default=1, min=1, max=5, step=1, precision=0, scale=1, unit="bins"),
                              group='Camera',
                              tooltip="Vertical pixel bin size to set on camera.\n"
                                      "Note: this is software binning, as opposed to hardware binning.")

        # relevant devices - sinara
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('scheduler')
        self.readout_subsequence = Readout(self, time_readout_us=3000)

        # relevant devices - labrad
        self.setattr_device('shutters')
        self.setattr_device('oven')
        self.setattr_device('aperture')
        self.setattr_device('trap_dc')
        self.setattr_device('camera')
        self.setattr_device('flipper')

        # define magic numbers
        # todo: document these
        self.pmt_dark_threshold_counts = 15
        self.time_runtime_max_s = 480.
        self.time_aramp_pulse_s = 2.
        self.path_image_save = r"\\eric.physics.ucla.edu\groups\motion\Data"

        self.pmt_sample_num = 100
        self.pmt_flip_threshold_std = 2
        self.cam_sample_num = 10
        self.cam_flip_threshold_avg = 0.15

    def prepare_experiment(self):
        """
        Prepare experimental values and precompute/preallocate to reduce kernel overheads.
        """
        '''
        HARDWARE VALUES
        '''
        self.start_time_s = 0.  # store runtime to check timeout

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


        '''
        FILEPATHS
        '''
        # todo: save image data to a dataset as well so we can image dynamically
        # todo: migrate data storage to a dataset so we can re-process images later
        year_month =        datetime.today().strftime('%Y-%m')
        year_month_day =    datetime.today().strftime('%Y-%m-%d')
        self.data_path = os.path.join(self.path_image_save, year_month, year_month_day)
        # ensure data folder exists; create it if it doesn't exist
        os.makedirs(self.data_path, exist_ok=True)


    '''
    INITIALIZATION & CLEANUP
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # set readout profile for beams and turn them on
        self.pump.readout()
        self.pump.on()
        self.repump_qubit.on()
        self.repump_cooling.on()

        # deterministically set flipper to camera and synchronize timeline
        self.set_flipper_to_camera()
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

    @rpc
    def initialize_labrad_devices(self) -> TNone:
        """
        Initialize labrad devices via RPC.
        """
        print("\tINITIALIZE - BEGIN")
        # note: try/except block must be here b/c initialize_experiment can't have
        #   try/except (b/c it's a kernel_from_string, which makes try/except difficult)
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
            print("Error during initialize_labrad_devices: {:}".format(repr(e)))
            self.cleanup_devices()

        finally:
            print("\tINITIALIZE - FINISH")

    @rpc
    def cleanup_devices(self) -> TNone:
        """
        Set all devices to states as if ion was loaded.
        All but trap electrodes set to original state --- trap electrodes set to final trapping potential.
        """
        print("\tCLEANUP - BEGIN")
        # turn off oven
        self.oven.set_oven_voltage(0.)
        self.oven.set_oven_current(0.)
        self.oven.toggle(False)

        # close shutters
        self.shutters.toggle_377_shutter(False)
        self.shutters.toggle_423_shutter(False)

        # set trap parameters for normal operation
        self.trap_dc.set_h_shim_voltage(self.end_h_shim_voltage)
        self.trap_dc.set_v_shim_voltage(self.end_v_shim_voltage)
        self.trap_dc.set_aramp_voltage(self.end_aramp_voltage)
        # turn on the other DC channels
        self.trap_dc.h_shim_toggle(True)
        self.trap_dc.v_shim_toggle(True)
        self.trap_dc.aramp_toggle(True)
        # ramp endcaps to normal values
        self.trap_dc.ramp_both_endcaps([self.end_east_endcap_voltage, self.end_west_endcap_voltage],
                                       [100., 100.])
        # note: we add sleep and set voltages AGAIN to reflect updates on GUIs
        sleep(2)
        self.trap_dc.set_east_endcap_voltage(self.end_east_endcap_voltage)
        self.trap_dc.set_west_endcap_voltage(self.end_west_endcap_voltage)

        # close aperture
        # note: do this last to keep ion happy during voltage adjustments
        self.aperture.close_aperture()

        if self.set_to_pmt_after_loading: self.flipper.flip()
        print("\tCLEANUP - FINISH")


    """
    MAIN LOGIC
    """
    @rpc
    def run_main(self) -> TNone:
        """
        Run till ion is loaded or timeout.
        """
        num_ions = 0    # run loading loop until we load desired_num_of_ions
        self.start_time_s = time()  # get start time so we know if we exceed max time

        try:
            while (num_ions != self.desired_num_of_ions) and (num_ions != -1):
                self.check_termination()    # periodically check termination

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

        # note: no need to cleanup if TerminationRequested b/c check_termination does it automatically
        except TerminationRequested:
            print("Termination Requested - terminating.")

        # catch all other errors
        except Exception as e:
            print("Error: {:}\nStopping execution & cleaning up.".format(repr(e)))
            self.cleanup_devices()

        # if no error, simply clean up
        else:   self.cleanup_devices()

    @rpc
    def load_ion(self) -> TInt32:
        """
        Function for loading an ion.
        This function runs continuously and blocks, stopping only when either sufficient ions are loaded,
            or until we time out.
        Returns:
            TInt32: number of detected ions in the trap
        """
        print("\tLOAD ION - BEGIN:")
        ion_spottings = 0   # improve detection certainty by lowpass filtering detection events

        # run loop while we don't see ions
        while True:
            # extract number of ions on camera
            num_ions = self.process_image('pre_aramp_original.png', 'pre_aramp_manipulated.png')

            # check if we've reached a stop condition (e.g. max_time reached, termination_requested)
            if (time() - self.start_time_s) > self.time_runtime_max_s:
                print("\tPROBLEM: TOOK OVER 15 MIN TO LOAD - ENDING PROGRAM")
                return -1
            else:   self.check_termination()

            # check to see if we have loaded enough ions
            if num_ions >= self.desired_num_of_ions:
                ion_spottings += 1
                # ensure camera sees ion in 3 consecutive images to prevent singular false positive
                if ion_spottings >= 3:
                    print("\t\t{:d} ION(s) LOADED".format(num_ions))
                    return num_ions
            else:   ion_spottings = 0  # reset if image analysis shows no ions in trap

        return 0

    @rpc
    def aramp_ions(self) -> TInt32:
        """
        Pulse A-ramp to get rid of excess ions.
        The A-ramp voltage progresses with aramp_ions_voltage_list until either the excess ions are ejected,
            or we reach the end of the list.
        This function runs continuously and blocks, stopping only when sufficient ions are jected,
            or until we time out.
        Returns:
            TInt32: Number of ions.
        """
        # todo: catch error conditions
        self.aperture.open_aperture()

        print("\tARAMP EJECTION - START")
        for aramp_voltage in self.aramp_ions_voltage_list:
            self.check_termination()

            # pulse A-ramp for given period
            print("\t\tARAMPING @ VOLTAGE: {:.1f} V".format(aramp_voltage))
            self.trap_dc.set_aramp_voltage(aramp_voltage)
            sleep(self.time_aramp_pulse_s)
            self.trap_dc.set_aramp_voltage(self.end_aramp_voltage)
            sleep(self.time_aramp_pulse_s)

            # get number of ions
            num_ions = self.process_image('post_aramp_original.png', 'post_aramp_manipulated.png')
            print("\t\tIONS REMAINING: {:d}".format(num_ions))
            if num_ions <= self.desired_num_of_ions:
                return num_ions

        self.aperture.close_aperture()
        return 0

    @rpc
    def process_image(self, filepath1: TStr="original.png", filepath2: TStr="manipulated.png") -> TInt32:
        """
        Process image data from camera and extract the number of ions in the trap.
        Arguments:
            filepath1   (string)    : save path for the original (no post-processing) image.
            filepath2   (string)    : save path for the processed image.
        Returns:
            TInt32: Number of ions in the trap
        """
        self.update_completion_monitor(time() - self.start_time_s)  # update completion monitor

        # get camera data and reshape into image
        image_arr = self.camera.get_most_recent_image()
        data = np.reshape(image_arr, (self.image_width_pixels, self.image_width_pixels))
        imsave(os.path.join(self.data_path, filepath1), data)

        # threshold & rescale data
        # todo: set 1000 as some parameter for min scatter value
        data *= data > 1000
        data = np.uint8(((data - np.min(data)) / (np.max(data) - np.min(data))) * 255)
        # use only upper 1% quantile of data
        data *= data > np.quantile(data, 0.99)
        imsave(os.path.join(self.data_path, filepath2), data)

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
        return len(unique_locs)


    """
    HELPER FUNCTIONS
    """
    @kernel(flags={"fast-math"})
    def set_flipper_to_camera(self) -> TNone:
        """
        Set the flipper to send light to the camera.
        """
        ### PREPARE 397nm (BEAM AND APERTURE TO CHARACTERIZE SCATTER ###
        # ensure aperture is open w/397nm off to ensure camera/PMT configuration is successful
        self.core.break_realtime()
        self.pump.beam.set_att(31.5 * dB)   # have to set att b/c readout.run() turns on beams
        self.aperture.open_aperture()
        sleep(2)    # add delay to ensure aperture closes successfully


        ### check PMT and camera status ###
        # get PMT and camera counts with aperture open
        self.core.wait_until_mu(now_mu())
        pmt_start_avg, pmt_start_std = self._burst_pmt()
        cam_start_avg, cam_start_std = self._burst_cam()
        print('\t\tpmt w/AP OPEN:', pmt_start_avg, '\t', pmt_start_std)
        print('\t\tcam w/AP OPEN:', cam_start_avg, '\t', cam_start_std)

        # get PMT and camera counts with aperture closed
        self.core.wait_until_mu(now_mu())
        self.aperture.close_aperture()
        sleep(2)    # add delay to ensure aperture closes successfully
        self.core.wait_until_mu(now_mu())
        pmt_stop_avg, pmt_stop_std = self._burst_pmt()
        cam_stop_avg, cam_stop_std = self._burst_cam()
        print('\t\tpmt w/397 ON:', pmt_stop_avg, '\t', pmt_stop_std)
        print('\t\tcam w/397 ON:', cam_stop_avg, '\t', cam_stop_std)

        # return 397nm att to normal
        self.core.break_realtime()
        self.pump.beam.set_att_mu(self.pump.att_pump_mu)


        ### DETERMINE FLIPPER STATUS AND FLIP ACCORDINGLY ###
        # determine whether flipper set to camera or PMT based on ability to register 397nm on/off
        #   above some threshold (determined as a multiple of the stdev)
        flip_status_pmt = (
                abs(pmt_start_avg - pmt_stop_avg) >
                self.pmt_flip_threshold_std * min(pmt_start_std, pmt_stop_std)
        )
        flip_status_cam = (
                abs(cam_start_avg - cam_stop_avg) >
                self.cam_flip_threshold_avg * min(cam_start_avg, cam_stop_avg)
        )
        # tmp remove
        print('\t\tval pmt:', abs(pmt_start_avg - pmt_stop_avg))
        print('\t\tthreshold pmt:', self.pmt_flip_threshold_std * min(pmt_start_std, pmt_stop_std))
        print('\t\tval cam:', abs(cam_start_avg - cam_stop_avg))
        print('\t\tthreshold cam:', self.cam_flip_threshold_avg * min(cam_start_avg, cam_stop_avg))
        print('\t\tflip status pmt:', flip_status_pmt)
        print('\t\tflip status cam:', flip_status_cam)
        # tmp remove

        # flip flipper according to status
        if flip_status_pmt and (not flip_status_cam):   # flipper set to PMT => flip to cam
            self.flipper.flip()
        elif (not flip_status_pmt) and flip_status_cam: # flipper already on cam => flip to PMT
            pass
        elif flip_status_pmt == flip_status_cam:        # error condition
            raise ValueError("PROBLEM: unable to determine whether flipper set to camera or PMT.")

    @kernel(flags={"fast-math"})
    def _burst_pmt(self) -> TTuple([TFloat, TFloat]):
        """
        Record PMT counts in a burst.
        :return: avg and std of the PMT counts (i.e. tuple(counts_avg, counts_std))
        """
        counts = [0] * self.pmt_sample_num  # store counts

        # store PMT count events
        self.core.break_realtime()
        for i in range(self.pmt_sample_num):
            self.readout_subsequence.run()
            delay_mu(128)
        self.core.break_realtime()

        # retrieve counts from PMT
        for i in range(self.pmt_sample_num):
            counts[i] = self.readout_subsequence.fetch_count()
            delay_mu(128)

        # calculate count statistics - mean
        counts_avg = 0.
        for val in counts:
            counts_avg += val
        counts_avg /= self.pmt_sample_num

        # calculate count statistics - stdev
        counts_std = 0.
        for val in counts:
            counts_std += (val - counts_avg) ** 2
        counts_std = np.sqrt(counts_std / (self.pmt_sample_num - 1))
        return counts_avg, counts_std

    @rpc
    def _burst_cam(self) -> TTuple([TFloat, TFloat]):
        """
        Record camera counts in a burst.
        :return: avg and std of the camera counts.
        """
        cam_vals = np.zeros(self.cam_sample_num)
        for i in range(self.cam_sample_num):
            cam_vals[i] = np.sum(self.camera.get_most_recent_image())

        return np.mean(cam_vals), np.std(cam_vals)

    @rpc
    def check_termination(self) -> TNone:
        """
        OVERRIDE base check_termination to ensure devices always clean up
        in case of termination.
        """
        if self.scheduler.check_termination():
            self.cleanup_devices()
            raise TerminationRequested

    @rpc(flags={"async"})
    def update_completion_monitor(self, elapsed_time_s: TFloat) -> TNone:
        """
        Update the completion monitor to inform user of elapsed loading time.
        :param elapsed_time_s: elapsed time in seconds.
        """
        self.set_dataset('management.dynamic.completion_pct',
                         round(elapsed_time_s / self.time_runtime_max_s * 100., 3),
                         broadcast=True, persist=True, archive=False)

