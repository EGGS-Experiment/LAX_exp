import numpy as np
from artiq.experiment import *
from matplotlib.pyplot import imsave

import os
from time import time, sleep
from datetime import datetime

from LAX_exp.language import *
from LAX_exp.system.subsequences import Readout

class CollisionDetection(LAXExperiment, Experiment):
    """
    Utility: Collision Detection

    Automatically tell if a collision has occurred and open aperture and lower endcaps to recrystallize chain
    """
    name = 'Collision Detection'

    kernel_invariants = {
        # hardware values

        # magic numbers
        "pmt_sample_num", "pmt_dark_threshold_counts", "pmt_flip_threshold_std",
        "cam_sample_num", "cam_flip_threshold_avg",

        # camera configuration
        "image_region", "path_image_save", "data_path",
    }

    def build_experiment(self):

        # timing arguments
        group_name_timing = 'Timing'
        self.setattr_argument('time_run_hrs',
                              NumberValue(default=8,
                                          step=0.01, precision=2,
                                          min = 0.01, max=12.,
                                          scale=1, unit = "hrs"),
                              group = group_name_timing)

        # recrystallization arguements
        group_name_recrystallization = 'Recrystallization'
        self.setattr_argument('max_recrystallization_attempts',
                                                                    NumberValue(default=10,precision=0,
                                                                                step=1, type='int',
                                                                                max = 100, min=1,
                                                                                ),
                                                                    group = group_name_recrystallization)
        self.setattr_argument('time_recrystallization_s',
                                                         NumberValue(default=10,precision=2,
                                                                     step = 0.1, max = 50.,
                                                                     min = 0.1, unit = 's'),
                                                         group=group_name_recrystallization)


        # low secular trap arguments
        group_name_low_sec_trap = 'Low Trap Secular Params'
        self.setattr_argument('low_east_endcap_voltage',
                              NumberValue(default=18., precision=1,
                                          step=0.1, min=0., max=400., scale=1., unit="V"),
                              group=group_name_low_sec_trap)
        self.setattr_argument('low_west_endcap_voltage',
                              NumberValue(default=20., precision=1, step=0.1,
                                          min=0., max=400., scale=1., unit="V"),
                              group=group_name_low_sec_trap)

        # high secular trap arguments
        group_name_high_sec_trap = 'High Trap Secular Params'

        self.setattr_argument('high_east_endcap_voltage',
                              NumberValue(default=218., precision=1, step=0.1,
                                          min=0., max=400., scale=1., unit="V"),
                              group=group_name_high_sec_trap)


        self.setattr_argument('high_west_endcap_voltage',
                              NumberValue(default=320., precision=1, step=0.1,
                                         min=0., max=400., scale=1., unit="V"),
                              group=group_name_high_sec_trap)

        self.setattr_argument('h_shim_voltage',
                              NumberValue(default=29.,
                            precision=1, step=0.1, min=0.,
                            max=150., scale=1., unit="V"),
                              group=group_name_high_sec_trap)

        self.setattr_argument('v_shim_voltage',
                              NumberValue(default=69, precision=1,
                                          step=0.1, min=0., max=150.,
                                          scale=1., unit="V"),
                              group=group_name_high_sec_trap)

        self.setattr_argument('aramp_voltage',
                              NumberValue(default=2.8, precision=1,
                                          step=0.1, min=0., max=50.,
                                          scale=1., unit="V"),
                              group=group_name_high_sec_trap)

        # 397 main pass (pump) parameters
        group_name_397 = '397_parameters'
        self.setattr_argument("freq_pump_mhz", NumberValue(100,
                                                           min=60, max=200, step=0.01,
                                                          unit="MHz", scale=1, precision=6
                                                    ),
                                                    group=group_name_397)

        self.setattr_argument("ampl_pump_pct",
                              NumberValue(default=18., precision=3, step=5,
                                          min=11, max=50, scale=1., unit='%'),
                              group=group_name_397)

        self.setattr_argument("att_pump_db",
                              NumberValue(default=14., precision=1, step=0.5, min=11.,
                                          max=31.5, scale=1., unit='dB'),
                              group=group_name_397)

        # imaging parameters
        group_name_cam = 'Camera'
        self.setattr_argument('time_exposure_s',
                                NumberValue(
                                    default=0.1, min = 0.1, max=0.5,
                                    step = 0.1, precision=1, scale=1, unit = 's'
                                ),
                              group = group_name_cam
                                )
        self.setattr_argument('image_center_x', NumberValue(default=256, min=1, max=512,
                                                            step=50, precision=0, scale=1, unit="pixels"),
                              group=group_name_cam,
                              tooltip="Image center position (x).")
        self.setattr_argument('image_center_y', NumberValue(default=256, min=1, max=512,
                                                            step=50, precision=0, scale=1, unit="pixels"),
                              group=group_name_cam,
                              tooltip="Image center position (y).")
        self.setattr_argument('image_width_pixels', NumberValue(default=350, min=100, max=450,
                                                                step=50, precision=0, scale=1, unit="pixels"),
                              group=group_name_cam,
                              tooltip="Square width of the total 512x512 image region on camera.\n"
                                      "Defined relative to (image_center_x, image_center_y).\n"
                                      "Note: too large an image width (~450 pixels) may contain scatter and confuse the analysis.")
        self.setattr_argument('horizontal_binning', NumberValue(default=1, min=1, max=5, step=1, precision=0, scale=1, unit="bins"),
                              group=group_name_cam,
                              tooltip="Horizontal pixel bin size to set on camera.\n"
                                      "Note: this is software binning, as opposed to hardware binning.")
        self.setattr_argument('vertical_binning',   NumberValue(default=1, min=1, max=5, step=1, precision=0, scale=1, unit="bins"),
                              group=group_name_cam,
                              tooltip="Vertical pixel bin size to set on camera.\n"
                                      "Note: this is software binning, as opposed to hardware binning.")

        # relevant devices - sinara
        self.setattr_device('pump')
        self.setattr_device('probe')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('scheduler')
        self.readout_subsequence = Readout(self, time_readout_us=3000)

        # relevant devices - labrad
        self.setattr_device('aperture')
        self.setattr_device('trap_dc')
        self.setattr_device('camera')
        self.setattr_device('flipper')

        # define magic numbers
        self.pmt_dark_threshold_counts = 5
        self.path_image_save = r"\\eric.physics.ucla.edu\groups\motion\Data"

        self.pmt_sample_num = 100
        self.pmt_flip_threshold_std = 2
        self.cam_sample_num = 10
        self.cam_flip_threshold_avg = 0.1

        # parameters for 397 collision event - leave hardcoded for now
        self.freq_collision_pump_mhz = 90
        self.att_collision_pump_db = 14
        self.ampl_collision_pump_pct = 18.

    def prepare_experiment(self):
        """
        Prepare experimental values and precompute/preallocate to reduce kernel overheads.
        """

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
        year_month =        datetime.today().strftime('%Y-%m')
        year_month_day =    datetime.today().strftime('%Y-%m-%d')
        self.data_path = os.path.join(self.path_image_save, year_month, year_month_day)
        self.recovery_videos_path = os.path.join(self.data_path, "recovery_videos")
        # ensure data folder exists; create it if it doesn't exist
        os.makedirs(self.data_path, exist_ok=True)
        os.makedirs(self.recovery_videos_path, exist_ok=True)

        self.rest_time_mu = self.core.seconds_to_mu(5)

        ### CONVERT LASER VALUES TO MACHINE UNITS ###
        self.freq_pump_ftw = self.pump.frequency_to_ftw(self.freq_pump_mhz*MHz)
        self.ampl_pump_asf = self.pump.amplitude_to_asf(self.ampl_pump_pct / 100.)
        self.att_pump_mu = att_to_mu(self.att_pump_db * dB)

        self.freq_collision_pump_ftw = self.pump.frequency_to_ftw(self.freq_collision_pump_mhz*MHz)
        self.ampl_collison_pump_asf = self.pump.amplitude_to_asf(self.ampl_collision_pump_pct / 100.)
        self.att_collision_pump_mu = att_to_mu(self.att_collision_pump_db * dB)

        self.time_recrystallization_mu = self.core.seconds_to_mu(self.time_recrystallization_s)

        # convert run time from hours to seconds and then machine units
        self.time_run_mu = self.core.seconds_to_mu(self.time_run_hrs * 3600)

        # timing requirements for rpc
        self.time_voltage_ramp_wait_s = 2.
        self.time_voltage_ramp_wait_mu = self.core.seconds_to_mu(self.time_voltage_ramp_wait_s)

        self.aperture_movement_wait_s = 2.
        self.aperture_movement_wait_mu = self.core.seconds_to_mu(self.aperture_movement_wait_s)

        self.succesful_recovery_attempts = 0
        self.recording_start_s = 0

    '''
    INITIALIZATION & CLEANUP
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # set readout profile for beams and turn them on
        self.pump.readout()
        self.pump.on()
        self.probe.off()
        self.repump_qubit.on()
        self.repump_cooling.on()

        # set pump to normal operating parameters
        self.pump.set_profile(self.pump.profile_readout)
        self.pump.set_att_mu(self.att_pump_mu)
        self.pump.set_mu(
            ftw = self.freq_pump_ftw,
            asf = self.ampl_pump_asf,
            profile = self.pump.profile_readout
        )

        self.initialize_labrad_devices()
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

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
            # Configure everything while acquisition is stopped.
            self.camera.stop_acquisition()
            self.camera.set_image_region(self.image_region)
            self.camera.set_exposure_time(self.time_exposure_s)
            self.camera.set_mode_shutter("Open")
            self.camera.set_mode_acquisition("Run till abort")
            self.camera.start_acquisition()

            '''SET UP TRAP'''
            # set endcaps to loading voltages
            self.trap_dc.set_east_endcap_voltage(self.high_east_endcap_voltage)
            self.trap_dc.set_west_endcap_voltage(self.high_west_endcap_voltage)

        except Exception as e:
            print("Error during initialize_labrad_devices: {:}".format(repr(e)))
            raise NotImplementedError('Issue initialized cleanup devices. Cleaning up.')

        finally:
            print("\tINITIALIZE - FINISH")

    @rpc
    def cleanup_devices(self) -> TNone:
        """
        Set all devices to states as if ion was loaded.
        All but trap electrodes set to original state --- trap electrodes set to final trapping potential.
        """
        print("\tCLEANUP - BEGIN")

        # set trap parameters for normal operation
        self.trap_dc.set_h_shim_voltage(self.h_shim_voltage)
        self.trap_dc.set_v_shim_voltage(self.v_shim_voltage)
        self.trap_dc.set_aramp_voltage(self.aramp_voltage)
        # turn on the other DC channels
        self.trap_dc.h_shim_toggle(True)
        self.trap_dc.v_shim_toggle(True)
        self.trap_dc.aramp_toggle(True)
        # ramp endcaps to normal values
        self.trap_dc.ramp_both_endcaps([self.high_east_endcap_voltage, self.high_west_endcap_voltage],
                                       [100., 100.])
        # note: we add sleep and set voltages AGAIN to reflect updates on GUIs
        sleep(2)
        self.trap_dc.set_east_endcap_voltage(self.high_east_endcap_voltage)
        self.trap_dc.set_west_endcap_voltage(self.high_west_endcap_voltage)

        # close aperture
        # note: do this last to keep ion happy during voltage adjustments
        self.aperture.close_aperture()

    """
    MAIN LOGIC
    """
    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        """
        Run till timeout or failed collision recrystallizion
        """

        # get initial background counts
        # determine counts with laser on
        self.pump.on()

        # ensure aperture is closed
        self.aperture.close_aperture()
        delay_mu(self.aperture_movement_wait_mu)
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

        print('\tGETTING BRIGHT IMAGES')
        cam_vals_bright = self._burst_cam()
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

        # ensure wall clock has caught up to scheduler to we dont schedule pump
        # turning on too far in future and taking cam pictures before it is turned on
        self.pump.off()
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

        print('\tGETTING BACKGROUND IMAGES')
        cam_vals_bg = self._burst_cam()
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

        print('\tFINISHED COLLECTING BRIGHT AND BACKGROUND IMAGES')

        cam_threshold = self.find_camera_threshold(
            cam_vals_bright,
            cam_vals_bg,
        )

        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

        time_start_mu = self.core.get_rtio_counter_mu()

        # turn pump back on
        self.pump.on()

        try:
            while (self.core.get_rtio_counter_mu() - time_start_mu) < self.time_run_mu:
                self.check_termination()    # periodically check termination
                self.core.break_realtime()

                # grab camera vals
                cam_vals = self._burst_cam()
                self.core.wait_until_mu(now_mu())
                self.core.break_realtime()
                mean_cam_vals = 0.0
                for cam_val in cam_vals:
                    mean_cam_vals += cam_val
                mean_cam_vals /= len(cam_vals)

                if mean_cam_vals < cam_threshold:
                    # detected collision
                    print('\tCOLLISION DETECTED. '
                          '\n\tCOLLISION NUMBER: ', self.succesful_recovery_attempts+1,
                          '\n\tCAMERA COUNTS: ', mean_cam_vals)
                    crystallization_success = self.crystallize_chain_after_collision(
                        camera_threshold=cam_threshold,
                        time_start_mu = time_start_mu
                    )
                    if crystallization_success:
                        print('\tIONS RECRYSTALLIZED SUCCESSFULLY')
                    if not crystallization_success:
                        raise TerminationRequested("Termination Requested - Unable to recrystallize ion chain")

                self.update_completion_monitor(self.core.mu_to_seconds(
                    self.core.get_rtio_counter_mu() - time_start_mu
                )
                )  # update completion monitor

        # note: no need to cleanup if TerminationRequested b/c cleanup_experiment does it automatically
        except TerminationRequested:
            print("Termination Requested - terminating.")

        # catch all other errors
        except Exception as e:
            print("Error: \nStopping execution & cleaning up.")
            raise


    """
    HELPER FUNCTIONS
    """

    @kernel(flags={"fast-math"})
    def crystallize_chain_after_collision(self,
                                          camera_threshold: TFloat,
                                          time_start_mu: TInt64) -> TBool:

        crystallization_success = False

        for attempt in range(self.max_recrystallization_attempts):

            """
            Follwing procedure to recrystallize chain

            1. set AOM to be more red-detuned than usual 
            2. open aperture to let red-detuned carrier beam
            3. lower endcaps to reduce endcaps 
            4. turn off shims and aramp voltages
            5. let sit for given time to recrystallize
            6. detect if chain has recrystallized via ion chain
            """

            self.check_termination()
            self.core.break_realtime()
            if (self.core.get_rtio_counter_mu() - time_start_mu) >= self.time_run_mu:
                break

            self.begin_recovery_recording()
            self.core.break_realtime()

            # set pump to collision parameters
            self.pump.set_mu(
                ftw = self.freq_collision_pump_ftw,
                asf = self.ampl_collison_pump_asf,
                profile=self.pump.profile_readout
            )
            self.pump.set_att_mu(self.att_collision_pump_mu)

            # open aperture
            self.aperture.open_aperture()
            delay_mu(self.aperture_movement_wait_mu)


            # ramp endcaps DOWN
            self.trap_dc.ramp_both_endcaps([self.low_east_endcap_voltage,
                                            self.low_west_endcap_voltage],
                                           [100., 100.])
            delay_mu(self.time_voltage_ramp_wait_mu)

            # turn off shim and aramp electrodes
            self.trap_dc.h_shim_toggle(False)
            self.trap_dc.v_shim_toggle(False)
            self.trap_dc.aramp_toggle(False)

            # let timeline catch up
            self.core.wait_until_mu(now_mu())
            self.core.break_realtime()

            # try to crystallize for provided time
            delay_mu(self.time_recrystallization_mu)

            # set pump to normal parameters
            self.pump.set_mu(
                ftw = self.freq_pump_ftw,
                asf = self.ampl_pump_asf,
                profile=self.pump.profile_readout
            )
            self.pump.set_att_mu(self.att_pump_mu)

            # close aperture and wait for timeline to catch up
            self.core.wait_until_mu(now_mu())
            self.core.break_realtime()
            self.aperture.close_aperture()
            delay_mu(self.aperture_movement_wait_mu)
            self.core.wait_until_mu(now_mu())
            self.core.break_realtime()

            # take camera vals to see if ions have recrystallized
            cam_vals = self._burst_cam()
            mean_cam_vals = 0.0
            for cam_val in cam_vals:
                mean_cam_vals += cam_val
            mean_cam_vals /= len(cam_vals)
            self.core.wait_until_mu(now_mu())
            self.core.break_realtime()

            if (mean_cam_vals > camera_threshold):
                crystallization_success = True

                # ramp endcaps back up and turns shims and aramp voltages back on
                self.trap_dc.ramp_both_endcaps([self.high_east_endcap_voltage,
                                                self.high_west_endcap_voltage],
                                               [100., 100.])
                delay_mu(self.time_voltage_ramp_wait_mu)
                self.trap_dc.h_shim_toggle(True)
                self.trap_dc.v_shim_toggle(True)
                self.trap_dc.aramp_toggle(True)
                self.core.wait_until_mu(now_mu())
                self.core.break_realtime()

                # take camera vals to see if ions are still crystallized after raising endcaps
                cam_vals = self._burst_cam()
                mean_cam_vals = 0.0
                for cam_val in cam_vals:
                    mean_cam_vals += cam_val
                mean_cam_vals /= len(cam_vals)
                self.core.wait_until_mu(now_mu())
                self.core.break_realtime()

                if mean_cam_vals > camera_threshold:
                    crystallization_success = True

                    self.finish_recovery_recording(
                        crystallization_success,
                    )
                    self.core.break_realtime()
                    return crystallization_success
                else:
                    crystallization_success = False

                self.finish_recovery_recording(crystallization_success)
                self.core.break_realtime()

        return crystallization_success

    @rpc
    def find_camera_threshold(self,
                              cam_vals_bright,
                              cam_vals_bg) -> TFloat:
        bright = float(np.mean(cam_vals_bright))
        background = float(np.mean(cam_vals_bg))

        # heuristic method - determine if there is a better way
        threshold = background + 0.35 * (bright - background)

        print("\tInitial bright:", bright)
        print("\tBackground:", background)
        print("\tThreshold:", threshold)

        return float(threshold)

        return threshold

    @rpc
    def determine_cam_stats(self,
                            cam_vals) -> TTuple([TFloat, TFloat]):

        return float(np.mean(cam_vals)), float(np.std(cam_vals))

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
        delay_mu(self.core.seconds_to_mu(3.))


        ### check PMT and camera status ###
        # get PMT and camera counts with aperture open
        self.core.wait_until_mu(now_mu())
        pmt_start_avg, pmt_start_std = self._burst_pmt()
        cam_start_vals = self._burst_cam()
        cam_start_avg, cam_start_std = self.determine_cam_stats(cam_start_vals)
        print('\t\tpmt w/AP OPEN:', pmt_start_avg, '\t', pmt_start_std)
        print('\t\tcam w/AP OPEN:', cam_start_avg, '\t', cam_start_std)

        # get PMT and camera counts with aperture closed
        self.core.wait_until_mu(now_mu())
        self.aperture.close_aperture()
        sleep(2)    # add delay to ensure aperture closes successfully
        self.core.wait_until_mu(now_mu())
        pmt_stop_avg, pmt_stop_std = self._burst_pmt()
        cam_stop_vals = self._burst_cam()
        cam_stop_avg, cam_stop_std = self.determine_cam_stats(cam_stop_vals)
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
    def _burst_cam(self) -> TArray(TFloat, 1):
        """
        Record camera counts in a burst.
        :return: cam vals
        """

        # discard transitional image
        self.check_termination()
        self.camera.wait_for_acquisition()
        self.camera.get_most_recent_image()

        self.camera.get_most_recent_image()
        cam_vals = np.zeros(self.cam_sample_num)
        for i in range(self.cam_sample_num):
            self.camera.wait_for_acquisition()
            image = np.asarray(
                self.camera.get_most_recent_image(),
                dtype=np.float64
            )

            background = np.percentile(image, 97)
            signal = np.maximum(image - background, 0.0)

            cam_vals[i] = np.sum(signal)

        return cam_vals

    @rpc
    def check_termination(self) -> TNone:
        """
        OVERRIDE base check_termination to ensure devices always clean up
        in case of termination.
        """
        if self.scheduler.check_termination():
            raise TerminationRequested

    @kernel(flags={"fast-math"})
    def cleanup_experiment(self) -> TNone:
        self.core.break_realtime()
        self.cleanup_devices()
        self.core.break_realtime()

    @rpc(flags={"async"})
    def update_completion_monitor(self, elapsed_time_s: TFloat) -> TNone:
        """
        Update the completion monitor to inform user of elapsed loading time.
        :param elapsed_time_s: elapsed time in seconds.
        """
        self.set_dataset('management.dynamic.completion_pct',
                         round(elapsed_time_s / (self.time_run_hrs * 3600) * 100., 3),
                         broadcast=True, persist=True, archive=False)


    """
    RPCs for video recording
    """

    @rpc
    def begin_recovery_recording(self) -> TNone:
        # Restart acquisition to begin a fresh recording.
        self.camera.stop_acquisition()
        # lengthen exposure time to allow longer video recording
        self.camera.set_exposure_time(2.5 * self.time_exposure_s)

        self.recording_start_s = time()
        self.recording_capacity = self.camera.camera.buffer_size()

        self.camera.start_acquisition()


    @rpc
    def finish_recovery_recording(self,
                                  successful: TBool
                                  ) -> TNone:
        self.camera.stop_acquisition()
        # return camera to regular exposure time
        self.camera.set_exposure_time(self.time_exposure_s)

        try:
            if not successful:
                return

            self.succesful_recovery_attempts += 1

            try:
                first, last = self.camera.buffer_new_images()
            except Exception as e:
                if "DRV_NO_NEW_DATA" in str(e):
                    print("No recovery images available.")
                    return
                raise

            frame_count = last - first + 1

            hbin, vbin, xstart, xend, ystart, yend = self.image_region
            width = (xend - xstart + 1) // hbin
            height = (yend - ystart + 1) // vbin

            frames = np.empty(
                (frame_count, height, width),
                dtype=np.uint32
            )

            batch_size = 20
            for start_frame in range(0, frame_count, batch_size):
                self.check_termination()

                num_images = min(batch_size, frame_count - start_frame)

                images = self.camera.get_oldest_images(num_images)
                frames[start_frame:start_frame + num_images] = np.asarray(
                    images, dtype=np.uint32
                ).reshape(num_images, height, width)

            filename = (
                f"{self.scheduler.rid}_recover_attempt_{self.succesful_recovery_attempts}_"
                f"{datetime.now():%Y_%m_%d_%H_%M_%S}.npz"
            )


            np.savez_compressed(
                os.path.join(self.recovery_videos_path, filename),
                frames=frames,
                image_region=np.asarray(self.image_region),
                recording_start_s=self.recording_start_s,
                first_frame_index=first,
                last_frame_index=last
            )

            print("Saved recovery recording:", filename)

        finally:
            # Resume acquisition for normal monitoring.
            self.camera.start_acquisition()
