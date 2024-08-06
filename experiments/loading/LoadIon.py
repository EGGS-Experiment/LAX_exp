import extensions
import numpy as np
from artiq.experiment import *
import skimage
from skimage import measure

from LAX_exp.base import LAXExperiment


class IonLoad(LAXExperiment, Experiment):
    """
    Experiment: Ion Load

    Gets the number of counts as a function of frequency for a fixed time.
    """
    name = 'Ion Load'

    def build_experiment(self):

        # laser arguments
        self.setattr_argument('freq_397_mhz',
                              NumberValue(default=110., ndecimals=1,
                                          step=0.1, min=90., max=120., unit="MHz"), group='397')
        self.setattr_argument('ampl_397', NumberValue(default=50., ndecimals=1,
                                                      step=0.1, min=0., max=100.), group='397')
        self.setattr_argument('att_397_dB', NumberValue(default=14., ndecimals=1,
                                                        step=0.1, min=0., max=31.5, unit="dB"), group='397')

        self.setattr_argument('freq_866_mhz', NumberValue(default=110., ndecimals=1,
                                                          step=0.1, min=90., max=120., unit="MHz"), group='866')
        self.setattr_argument('ampl_866', NumberValue(default=50., ndecimals=1,
                                                      step=0.1, min=0., max=100.), group='866')
        self.setattr_argument('att_866_dB', NumberValue(default=14., ndecimals=1,
                                                        step=0.1, min=0., max=31.5, unit="dB"), group='866')

        self.setattr_argument('freq_854_mhz',
                              NumberValue(default=110., ndecimals=1,
                                          step=0.1, min=90., max=120., unit="MHz"), group='854')
        self.setattr_argument('ampl_854', NumberValue(default=50., ndecimals=1,
                                                      step=0.1, min=0., max=100.), group='854')
        self.setattr_argument('att_854_dB', NumberValue(default=14., ndecimals=1,
                                                        step=0.1, min=0., max=31.5, unit="dB"), group='854')

        # pmt arguments
        self.setattr_argument('ion_count_threshold', NumberValue(default=120, ndecimals=0,
                                                                 step=1, min=80, max=250), group='Photon Counting')
        self.setattr_argument('pmt_sample_time_us',
                              NumberValue(default=3e-3, ndecimals=0,
                                          step=1, min=1e-6, max=5e-3, unit='us'), group='Photon Counting')

        # starting trap arguments
        self.setattr_argument('starting_east_endcap_voltage',
                              NumberValue(default=19., ndecimals=1, step=0.1, min=0., max=300.),
                              group='Starting Trap Parameters')
        self.setattr_argument('starting_west_endcap_voltage',
                              NumberValue(default=25., ndecimals=1, step=0.1, min=0., max=300.),
                              group='Starting Trap Parameters')

        # starting trap arguments
        self.setattr_argument('ending_east_endcap_voltage',
                              NumberValue(default=202., ndecimals=1, step=0.1, min=0., max=300.),
                              group='Ending Trap Parameters')
        self.setattr_argument('ending_west_endcap_voltage',
                              NumberValue(default=289., ndecimals=1, step=0.1, min=0., max=300.),
                              group='Ending Trap Parameters')
        self.setattr_argument('ending_v_shim_voltage',
                              NumberValue(default=66.9, ndecimals=1, step=0.1, min=0., max=150.),
                              group='Ending Trap Parameters')
        self.setattr_argument('ending_h_shim_voltage',
                              NumberValue(default=50.1, ndecimals=1, step=0.1, min=0., max=150.),
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

        # relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('shutters')
        self.setattr_device('gpp3060')
        self.setattr_device('pmt')
        self.setattr_device('aperture')
        self.setattr_device('trap_dc')
        self.setattr_device('scheduler')
        self.setattr_device('camera')

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

        # convert 397 parameters
        self.ftw_397 = extensions.mhz_to_ftw(self.freq_397_mhz)
        self.asf_397 = extensions.pct_to_asf(self.ampl_397)
        self.att_397 = extensions.att_to_mu(self.att_397_dB)

        # convert 854 parameters
        self.ftw_854 = extensions.mhz_to_ftw(self.freq_854_mhz)
        self.asf_854 = extensions.pct_to_asf(self.ampl_854)
        self.att_854 = extensions.att_to_mu(self.att_854_dB)

        # convert 866 parameters
        self.ftw_866 = extensions.mhz_to_ftw(self.freq_866_mhz)
        self.asf_866 = extensions.pct_to_asf(self.ampl_866)
        self.att_866 = extensions.att_to_mu(self.att_866_dB)

    @property
    def results_shape(self):
        return (2, 2)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):

        self.core.break_realtime()
        self.pump.beam.cpld.get_att_mu()

        # set 397 parameters
        self.pump.set_mu(self.ftw_397, asf=self.asf_397, profile=6)
        self.pump.set_att_mu(self.att_397)

        # set 854 parameters
        self.repump_qubit.set_mu(self.ftw_854, asf=self.asf_854, profile=6)
        self.repump_qubit.set_att_mu(self.att_854)

        # set 866 parameters
        self.repump_cooling.set_mu(self.ftw_866, asf=self.asf_866, profile=6)
        self.repump_cooling.set_att_mu(self.att_866)

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

        # open shutters
        self.shutters.open_377_shutter()
        self.shutters.open_423_shutter()
        self.core.break_realtime()

        # turn on the oven
        self.gpp3060.turn_oven_on()
        self.core.break_realtime()

        # open aperture
        self.aperture.open_aperture()
        self.core.break_realtime()

        # set camera region of interest
        self.camera.set_image_region(self.image_region)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        """
        Run till ion is loaded or timeout
        """
        self.core.break_realtime()  # add slack
        self.start_time = self.get_rtio_counter_mu()

        idx = 0
        breaker = False
        count_successes = 0
        num_ions = 0

        while count_successes < 10 or num_ions == 0:
            # increment loop counter
            idx += 1

            # read from pmt
            self.core.break_realtime()  # add slack
            self.pmt.count(self.pmt_sample_time_us)  # set pmt sample time
            counts = self.pmt.fetch_count()  # grab counts from PMT
            self.core.break_realtime()

            # read from camera
            self.camera.acquire_single_image()
            image = self.camera.get_most_recent_image()
            num_ions = self.get_num_ions(image)

            if num_ions > 0:
                self.print_ion_loaded_message(num_ions)


            if counts >= self.ion_count_threshold:
                count_successes += 1
                if count_successes >=10:
                    self.print_ion_loaded_message()

            if idx >= 49:
                idx = 0
                breaker = self.check_time()
                if breaker: break
                self.check_termination()  # check if termination is over threshold
                self.core.break_realtime()


    # ANALYSIS
    def analyze_experiment(self):

        self.cleanup_devices()

    @rpc
    def print_ion_loaded_message(num_ions = None):

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
        self.gpp3060.turn_oven_off()

        # close shutters
        self.shutters.open_377_shutter()
        self.shutters.open_423_shutter()

        # disconnect from labjack
        self.shutters.close_labjack()

        # close aperture
        self.aperture.close_aperture()

        # set trap parameters as if ion was loaded
        self.trap_dc.set_east_endcap_voltage(self.starting_east_endcap_voltage)
        self.trap_dc.set_west_endcap_voltage(self.starting_west_endcap_voltage)
        self.trap_dc.set_h_shim_voltage(self.ending_h_shim_voltage)
        self.trap_dc.set_v_shim_voltage(self.ending_v_shim_voltage)
        self.trap_dc.set_a_ramp2_voltage(self.ending_a_ramp2_voltage)

        # turn on the endcap channels
        self.trap_dc.east_endcap_on()
        self.trap_dc.west_endcap_on()
        self.trap_dc.h_shim_on()
        self.trap_dc.v_shim_on()
        self.trap_dc.a_ramp2_on()

        # ramp endcaps to end values
        self.trap_dc.ramp_both_endcaps([self.ending_east_endcap_voltage, self.ending_west_endcap_voltage],
                                       [100, 100])


    @rpc
    def get_num_ions(self, data) -> TInt32:
        """
        Return number of ions determined from camera image
        """

        # reshape 1D array into 2D array
        data = np.reshape(data, (self.image_width_pixels, self.image_height_pixels))

        # threshold data into binary image
        upper_percentile = np.percentile(data, 99)
        data[data < upper_percentile] = 0
        data[data >= upper_percentile] = 1

        # erode then dilate image to remove small bright regions (scatter) and accentuate larger regions (ions)
        kernel = np.ones((2, 2), np.uint8)
        for i in range(3):
            data = np.uint8(skimage.morphology.binary_erosion(data, kernel))
        for i in range(3):
            data = np.uint8(skimage.morphology.binary_dilation(data, kernel))

        # ensure binary image
        data[data > 0] = 1

        # label all pixels as belonging to a localized patter
        labels = measure.label(data)

        # return number of localized patter (-1 as background is given a label)
        return len(np.unique(labels)) - 1


    @kernel
    def check_time(self) -> TBool:
        """
        Check wall clock time to see if 10 min (600 sec) has elapsed with no ion loaded
        """
        self.core.break_realtime()
        return 600 > self.core.mu_to_seconds(
            self.get_rtio_counter_mu() - self.start_time)  # check if longer than 10 min

    @rpc
    def check_termination(self):
        """
        Check whether termination of the experiment has been requested.
        """
        if self.scheduler.check_termination():
            self.cleanup_devices()
            raise TerminationRequested
