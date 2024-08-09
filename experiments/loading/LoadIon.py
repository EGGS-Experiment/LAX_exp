import numpy as np
from artiq.experiment import *
import skimage

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from matplotlib import pyplot as plt


class IonLoad(LAXExperiment, Experiment):
    """
    Experiment: Ion Load

    Gets the number of counts as a function of frequency for a fixed time.
    """
    name = 'Ion Load'

    IMAGE_HEIGHT = 512
    IMAGE_WIDTH = 512
    PMT_SAMPLE_TIME_US = 3000
    PMT_SAMPLE_NUM = 30

    def build_experiment(self):

        # laser arguments
        self.setattr_argument('freq_397_mhz',
                              NumberValue(default=100., ndecimals=1,
                                          step=0.1, min=90., max=120., unit="MHz"), group='397')
        self.setattr_argument('ampl_397', NumberValue(default=50., ndecimals=1,
                                                      step=0.1, min=0., max=100.), group='397')
        self.setattr_argument('att_397_dB', NumberValue(default=14., ndecimals=1,
                                                        step=0.1, min=0., max=31.5, unit="dB"), group='397')

        self.setattr_argument('freq_866_mhz', NumberValue(default=112., ndecimals=1,
                                                          step=0.1, min=90., max=120., unit="MHz"), group='866')
        self.setattr_argument('ampl_866', NumberValue(default=23., ndecimals=1,
                                                      step=0.1, min=0., max=100.), group='866')
        self.setattr_argument('att_866_dB', NumberValue(default=14., ndecimals=1,
                                                        step=0.1, min=0., max=31.5, unit="dB"), group='866')

        self.setattr_argument('freq_854_mhz',
                              NumberValue(default=112., ndecimals=1,
                                          step=0.1, min=90., max=120., unit="MHz"), group='854')
        self.setattr_argument('ampl_854', NumberValue(default=18., ndecimals=1,
                                                      step=0.1, min=0., max=100.), group='854')
        self.setattr_argument('att_854_dB', NumberValue(default=14., ndecimals=1,
                                                        step=0.1, min=0., max=31.5, unit="dB"), group='854')

        # starting trap arguments
        self.setattr_argument('starting_east_endcap_voltage',
                              NumberValue(default=15., ndecimals=1, step=0.1, min=0., max=300.),
                              group='Starting Trap Parameters')
        self.setattr_argument('starting_west_endcap_voltage',
                              NumberValue(default=24., ndecimals=1, step=0.1, min=0., max=300.),
                              group='Starting Trap Parameters')

        # starting trap arguments
        self.setattr_argument('ending_east_endcap_voltage',
                              NumberValue(default=202., ndecimals=1, step=0.1, min=0., max=300.),
                              group='Ending Trap Parameters')
        self.setattr_argument('ending_west_endcap_voltage',
                              NumberValue(default=289., ndecimals=1, step=0.1, min=0., max=300.),
                              group='Ending Trap Parameters')
        self.setattr_argument('ending_v_shim_voltage',
                              NumberValue(default=69.2, ndecimals=1, step=0.1, min=0., max=150.),
                              group='Ending Trap Parameters')
        self.setattr_argument('ending_h_shim_voltage',
                              NumberValue(default=50.5, ndecimals=1, step=0.1, min=0., max=150.),
                              group='Ending Trap Parameters')
        self.setattr_argument('ending_a_ramp2_voltage',
                              NumberValue(default=2.0, ndecimals=1, step=0.1, min=0., max=100.),
                              group='Ending Trap Parameters')

        # image region parameters: MAX (450,450) TO PREVENT LASER SCATTER OFF ELECTRODES FROM CONFUSING ANALYSIS
        self.setattr_argument('image_width_pixels', NumberValue(default=400,
                                                                min=100, max=450, step=1,
                                                                scale=1, ndecimals=0), group='Camera')

        self.setattr_argument('image_height_pixels', NumberValue(default=400,
                                                                 min=100, max=450, step=1,
                                                                 scale=1, ndecimals=0), group='Camera')

        self.setattr_argument('horizontal_binning', NumberValue(default=1,
                                                                min=1, max=5, step=1,
                                                                scale=1, ndecimals=0), group='Camera')

        self.setattr_argument('vertical_binning', NumberValue(default=1,
                                                              min=1, max=5, step=1,
                                                              scale=1, ndecimals=0), group='Camera')

        self.setattr_argument('exposure_time_ms', NumberValue(default=100, ndecimals=1, min=1, max=1000,
                                                              step=0.1, unit='ms'), group='Camera')
        # relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        # self.setattr_device('shutters')
        self.setattr_device('oven')
        self.setattr_device('pmt')
        self.setattr_device('aperture')
        self.setattr_device('trap_dc')

        self.setattr_device('scheduler')
        self.setattr_device('camera')
        self.setattr_device('flipper')

    def prepare_experiment(self):

        # set up image region
        border_x = (self.IMAGE_WIDTH - self.image_width_pixels) / 2
        border_y = (self.IMAGE_HEIGHT - self.image_height_pixels) / 2

        start_x = int(border_x)
        end_x = int(self.IMAGE_WIDTH - border_x) - 1

        start_y = int(border_y)
        end_y = int(self.IMAGE_WIDTH - border_y) - 1

        self.image_region = (self.horizontal_binning, self.vertical_binning,
                             start_x, end_x, start_y, end_y)

        self.exposure_time_s = self.exposure_time_ms / ms

        # convert 397 parameters
        self.ftw_397 = self.pump.frequency_to_ftw(self.freq_397_mhz*MHz)
        self.asf_397 = pct_to_asf(self.ampl_397)
        self.att_397 = att_to_mu(self.att_397_dB)

        # convert 854 parameters
        self.ftw_854 = self.repump_qubit.frequency_to_ftw(self.freq_854_mhz*MHz)
        self.asf_854 = pct_to_asf(self.ampl_854)
        self.att_854 = att_to_mu(self.att_854_dB)

        # convert 866 parameters
        self.ftw_866 = self.repump_cooling.frequency_to_ftw(self.freq_866_mhz*MHz)
        self.asf_866 = pct_to_asf(self.ampl_866)
        self.att_866 = att_to_mu(self.att_866_dB)

        self.pmt_sample_time_mu = us_to_mu(self.PMT_SAMPLE_TIME_US)

        # initialize start time
        self.start_time = np.int64(0)

        self.set_dataset("counts", [])
        self.pmt_inverse_num_samples = 1 / self.PMT_SAMPLE_NUM

    @property
    def results_shape(self):
        return (2, 2)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):

        self.core.break_realtime()
        self.pump.beam.cpld.get_att_mu()
        self.core.break_realtime()

        # set 397 parameters
        self.pump.set_mu(self.ftw_397, asf=self.asf_397, profile=6)
        self.pump.set_att_mu(self.att_397)
        self.core.break_realtime()

        # set 854 parameters
        self.repump_qubit.set_mu(self.ftw_854, asf=self.asf_854, profile=6)
        self.repump_qubit.set_att_mu(self.att_854)
        self.core.break_realtime()

        # set 866 parameters
        self.repump_cooling.set_mu(self.ftw_866, asf=self.asf_866, profile=6)
        self.repump_cooling.set_att_mu(self.att_866)
        self.core.break_realtime()

        # turn on lasers
        self.pump.on()
        self.repump_qubit.on()
        self.repump_cooling.on()
        self.core.break_realtime()

        # set endcaps to loading voltages
        self.trap_dc.set_east_endcap_voltage(self.starting_east_endcap_voltage)
        self.trap_dc.set_west_endcap_voltage(self.starting_west_endcap_voltage)

        # turn on endcap channels and ensure others are off
        self.trap_dc.east_endcap_on()
        self.trap_dc.west_endcap_on()
        self.trap_dc.h_shim_off()
        self.trap_dc.v_shim_off()
        self.trap_dc.a_ramp2_off()
        self.core.break_realtime()

        # # open shutters
        # self.shutters.open_377_shutter()
        # self.shutters.open_423_shutter()
        # self.core.break_realtime()

        # turn on the oven
        self.oven.set_oven_voltage(1)
        self.oven.on()
        self.core.break_realtime()

        # open aperture
        self.aperture.open_aperture()
        self.core.break_realtime()

        # set camera region of interest and exposure time
        self.core.break_realtime()
        self.camera.set_image_region(self.image_region)
        self.core.break_realtime()

        self.set_flipper_to_pmt()
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        """
        Run till ion is loaded or timeout
        """

        self.core.break_realtime()  # add slack
        self.start_time = self.core.get_rtio_counter_mu()
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

        # define some variables for later use
        idx = 0
        breaker = False
        num_ions = 0
        camera_success = False
        ion_spottings = 0

        # gather counts from pmt to later verify flipper worked correctly
        delay_mu(100000)
        counts = self.sample_pmt_counts()
        self.append_to_dataset('counts', counts)
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

        # flip to camera
        self.flipper.flip()
        self.core.break_realtime()

        while camera_success is not True: # run loop while we don't see ions

            idx += 1

            # read from camera and gather number of ions
            self.camera.acquire_single_image()
            image = self.camera.get_most_recent_image()
            num_ions = self.show_ions(image)

            # every 5 image check if termination was requested or if we've taken too long to load
            if idx >= 5:
                self.check_termination()  # check if termination is over threshold
                self.core.break_realtime()
                idx = 0
                breaker = self.check_time()
                if breaker:
                    print("TOOK OVER 10 MIN TO LOAD --- ENDING PROGRAM")
                    break

            if num_ions > 0:
                ion_spottings += 1
                if ion_spottings >= 3: # ensure camera sees ion from 3 consecutive pictures to prevent singular false positive
                    camera_success = True
                    self.print_ion_loaded_message(num_ions)

            else:
                ion_spottings = 0 # reset if image analysis shows no ions in trap

            with parallel:
                self.check_termination()
                self.core.break_realtime()


    # ANALYSIS
    def analyze_experiment(self):
        print(self.get_dataset("counts"))
        self.cleanup_devices()

    @rpc
    def print_ion_loaded_message(self, num_ions=None):

        if num_ions is not None:
            print(f"{num_ions} IONS LOADED!!!")
        else:
            print("ION LOADED!!!")

    @rpc
    def cleanup_devices(self):
        """
        Set all devices to states as if ion was loaded
        All but trap electrodes set to original state
        """
        # turn off oven
        self.oven.set_oven_voltage(0)
        self.oven.off()

        # # close shutters
        # self.shutters.open_377_shutter()
        # self.shutters.open_423_shutter()

        # disconnect from labjack
        # self.shutters.close_labjack()

        # close aperture
        self.aperture.close_aperture()

        # set trap parameters as if ion was loaded
        self.trap_dc.set_east_endcap_voltage(self.starting_east_endcap_voltage)
        self.trap_dc.set_west_endcap_voltage(self.starting_west_endcap_voltage)
        self.trap_dc.set_h_shim_voltage(self.ending_h_shim_voltage)
        self.trap_dc.set_v_shim_voltage(self.ending_v_shim_voltage)
        self.trap_dc.set_a_ramp2_voltage(self.ending_a_ramp2_voltage)

        # turn on the endcap channels
        self.trap_dc.h_shim_on()
        self.trap_dc.v_shim_on()
        self.trap_dc.a_ramp2_on()

        # ramp endcaps to end values
        self.trap_dc.ramp_both_endcaps([self.ending_east_endcap_voltage, self.ending_west_endcap_voltage],
                                       [100., 100.])

        self.trap_dc.set_east_endcap_voltage(self.ending_east_endcap_voltage)
        self.trap_dc.set_west_endcap_voltage(self.ending_west_endcap_voltage)

    @kernel
    def check_time(self) -> TBool:
        """
        Check wall clock time to see if 10 min (600 sec) has elapsed with no ion loaded
        """
        self.core.break_realtime()
        return 600 < self.core.mu_to_seconds(
            self.core.get_rtio_counter_mu() - self.start_time)  # check if longer than 10 min


    @kernel
    def sample_pmt_counts(self):
        counts = 0.
        for i in range(self.PMT_SAMPLE_NUM):
            self.core.break_realtime()  # add slack
            self.pmt.count(self.pmt_sample_time_mu)  # set pmt sample time
            delay_mu(8)
            counts += self.pmt.fetch_count()  # grab counts from PMT
            self.core.break_realtime()

        counts = counts * self.pmt_inverse_num_samples
        delay_mu(100)
        return counts

    @kernel
    def set_flipper_to_pmt(self):
        self.flipper.flip()
        self.core.break_realtime()
        counts = self.sample_pmt_counts()

        if counts < 15:
            self.flipper.flip()
            self.core.break_realtime()
    @kernel
    def set_flipper_to_camera(self):
        self.flipper.flip()
        self.core.break_realtime()
        counts = 0.
        for i in range(self.pmt_num_samples):
            self.core.break_realtime()  # add slack
            self.pmt.count(self.pmt_sample_time_mu)  # set pmt sample time
            delay_mu(8)
            counts += self.pmt.fetch_count()  # grab counts from PMT
            self.core.break_realtime()

        counts = counts * self.pmt_inverse_num_samples
        delay_mu(100)  # see if this is enough

        if counts > 15:
            self.flipper.flip()
            self.core.break_realtime()

    @rpc
    def show_ions(self, data) -> TInt32:

        data = np.reshape(data, (self.image_width_pixels, self.image_height_pixels))

        plt.figure(1)
        plt.title("Original Image")
        plt.imshow(data)
        plt.savefig("Z:\motion\Pictures\original.png")

        upper_percentile = np.percentile(data, 99)
        data[data < upper_percentile] = 0
        data[data >= upper_percentile] = 1

        kernel = np.ones((2, 2), np.uint8)
        for i in range(3):
            data = np.uint8(skimage.morphology.binary_erosion(data, kernel))

        for i in range(3):
            data = np.uint8(skimage.morphology.binary_dilation(data, kernel))

        data[data > 0] = 1
        labels = skimage.measure.label(data)
        plt.figure(2)
        plt.imshow(data)
        plt.title("Manipulated Image")
        plt.savefig("Z:\motion\Pictures\manipulated.png")
        return len(np.unique(labels)) - 1
