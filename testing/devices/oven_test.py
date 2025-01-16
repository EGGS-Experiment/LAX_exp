
from artiq.experiment import *

import labrad
import numpy as np
from artiq.experiment import *

from LAX_exp.base import LAXExperiment

>>>>>>> 524534b4c277a6647269fd41893ccfcf493fad2b
class OvenTest(LAXExperiment, Experiment):
    """
    todo: document
    """
    name = "Oven Test"

    def build_experiment(self):

<<<<<<< HEAD
        # grab oven
=======
        # grab tickle device
>>>>>>> 524534b4c277a6647269fd41893ccfcf493fad2b
        self.setattr_device('oven')


    def prepare_experiment(self):
<<<<<<< HEAD

=======
>>>>>>> 524534b4c277a6647269fd41893ccfcf493fad2b
        pass

    @property
    def results_shape(self):
        return (2, 2)

<<<<<<< HEAD
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # ensure DMA sequences use profile 0
        self.dds_tickle.set_profile(0)
        self.dds_tickle.set_att_mu(self.att_tickle_mu)
        self.core.break_realtime()

    def run_main(self):

        # check if we still have HCL via laser scan
        laser_scan_sec_freq = 0.  # set up variable to take in secular frequency from laser scan
        chem_reaction = True  # assume chem reaction
        self._prepare_laser_scan_expid()  # prepare the arguments for the laser scan experiment
        rid_dj = self.scheduler.submit(pipeline_name='calibrations.Laser_Scan', expid=self.laser_scan_expid)
        while rid_dj in self.scheduler.get_status().keys():
            pass # hold until experiment is done running
        results = self.get_dataset("temp.laserscan.results")  # get the results of the laser scan and process
        peak_vals, _ = process_laser_scan_results(results)
        for i in range(100):  # run until we have HCl and not H2CL
            for peak_freq, peak_prob in peak_vals:
                if (self.expected_peak_freq_mhz - self.tolerable_laser_scan_drift_khz <= peak_freq
                        <= self.expected_lawhser_scan_peak_mhz + self.tolerable_laser_scan_drift_khz):
                    laser_scan_sec_freq_mhz = peak_freq * MHz
                    chem_reaction = False

            if chem_reaction:
                if i == 99:
                    raise TerminationRequested("Could Not Dissassociate Ion")
                self.pulse_aperture()
                results = self.submit_laser_scan()
                peak_vals, _ = process_laser_scan_results(results)

            else:
                break

        # check position of ion

        # get list of scanning frequencies and how long to chirp
        tickle_freqs_ftw_list = mhz_to_ftw(RangeScan(laser_scan_sec_freq - self.freq_tickle_span_kh/2 * kHz,
                                                    laser_scan_sec_freq+ self.freq_tickle_span_khz/2 * kHz,
                                                     self.tickle_chirp_steps))

        chirp_time_mu = self.core.seconds_to_mu((self.time_tickle_us * us/ self.tickle_chirp_steps))

        self.flip_mirror()  # flip mirror to ANDOR
        image, pixels_x, pixels_y = self.check_ion_position()  # grab bright regions
        xpos = self.analyze_image(image, pixels_x, pixels_y)  # find x pos of ion
        for i in range(100):
            if not (self.expected_x_pos - 1 <= xpos <= self.expected_x_pos + 1):
                # todo: is 397 on and should it be??? - damping force but would help keep ions???
                self.tickle_ion(tickle_freqs_ftw_list, chirp_time_mu) # tickle the ion with a chirped pulse
                image, pixels_x, pixels_y = self.check_ion_position()  # check if ion has moved back
                xpos = self.analyze_image(image, pixels_x, pixels_y)

            if i == 99 and not (self.expected_x_pos - 1 <= xpos <= self.expected_x_pos + 1):
                raise TerminationRequested("Could Not Tickle Ion Back to Proper Position")

        self.flip_mirror()  # flip mirror back to PMT

    def analyze(self):
        pass

    @kernel(flags={"fast-math"})
    def tickle_ion(self, tickle_freqs_ftw, chirp_time_mu):
        self.core.break_realtime()
        # todo: is 397 on and should it be???
        ampl_tickle_asf = np.int32(self.ampl_tickle_asf)

        for tickle_freq_ftw in tickle_freqs_ftw:
            self.dds_tickle.set_mu(tickle_freq_ftw, asf=ampl_tickle_asf, profile=0)  # configure tickle
            self.dds_tickle.on()  # tickle
            delay_mu(chirp_time_mu) # wait till time to move onto next freq

        self.dds_tickle.off()
        self.core.break_realtime()

        # support graceful termination
        with parallel:
            self.check_termination()
            self.core.break_realtime()

    def _prepare_laser_scan_expid(self):
        self.laser_scan_expid = {
                    "file":         "LAX_exp\\experiments\\LaserScan.py",
                    "class_name":   "LaserScan",
                    "log_level":    30,
                    "arguments": {
                        "repetitions":  self.repetitions,
                        "att_qubit_db": self.att_qubit_db,
                        "freq_qubit_scan_mhz": {
                            "center":       self.freq_check_ion_secular_freq_MHz,
                            "span":         self.freq_span_MHz,
                            "step":         self.freq_steo_MHz,
                            "randomize":    True,
                            "seed":         None,
                            "ty":           "CenterScan"
                        }
                    }
                }

    def check_ion_position(self):

        # identify_exposure_time = WithUnit(0.2, 's')
        # todo: figure out these parameters
        start_x = 1
        stop_x = 512
        start_y = 1
        stop_y = 512
        h_bin = stop_y - start_y + 1 # image becomes 1d array of just x direction
        v_bin = 5
        image_region = (h_bin,v_bin,start_x, stop_x, start_y, stop_y)

        pixels_x = (stop_x-start_x+1)/hbin
        pixels_y = (stop_y-start_y+1)/vbin

        # get single image from ANDOR
        self.cam.accquisition_stop()
        initial_exposure_time = self.cam.setup_exposure_time()
        # self.cam.setup_exposure_time(identify_exposure_time)
        initial_image_region = self.cam.get_image_region()
        self.cam.set_image_region(*image_region)
        self.cam.mode_shutter("Open")
        self.cam.mode_accquisition("Single Scan")
        self.cam.accquisition_start()
        self.cam.accquisition_wait()
        image = self.cam.accquire_image_recent()

        # return ANDOR to previous settings
        self.cam.setup_exposure_time(initial_exposure_time)
        self.cam.set_image_region(*initial_image_region)
        self.cam.mode_accquisition("Run till abort")
        self.cam.accquisition_start()
        self.cam.accquisition_wait()

        return image, pixels_x, pixels_y

    @kernel(flags={"fast-math"})
    def flip_mirror(self):
        self.ttl15.pulse(self.time_flipper_trigger_mu)
        delay_mu(self.time_flipper_wait_mu)

    def analyze_image(self, image, pixels_x, pixels_y):
        # todo: check if this accurately finds the x_pos of the center of the ion
        image = np.ravel(np.reshape(image, (pixels_x, pixels_y)))
        max_pos_x = np.argmax(image)
        return max_pos_x

    def pulse_aperture(self):
        self.aperture.open()
        self.aperture.wait(self.aperture_wait_time)
        self.aperture.close()
=======
    def initialize_experiment(self):
        pass

    def run_main(self):
        self.oven.set_oven_voltage(0.)
        self.oven.set_oven_current(0.)
        self.oven.toggle(False)
        print(self.oven.get_oven_current())
        print(self.oven.get_oven_voltage())

    def analyze_experiment(self):
        pass

>>>>>>> 524534b4c277a6647269fd41893ccfcf493fad2b

