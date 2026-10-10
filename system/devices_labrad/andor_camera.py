from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad
from os import environ


class AndorCamera(LAXDevice):
    """
    High-level API functions for using the Andor Camera
    """
    name = "camera"
    kernel_invariants = {
        "cxn", "camera"
    }

    def prepare_device(self):
        self.camera = self.cxn.andor_server

        # retrieve camera info - tuple(dimensions_x, dimensions_y)

    @rpc
    def acquire_detector_dimensions(self):
        return self.camera.info_detector_dimensions()

    @rpc
    def acquire_single_image(self, image_region=None, identify_exposure_time: TFloat=None) -> TNone:
        """
        Acquire a single image from the camera and then reset to previous setting
        :param image_region: area of ion trap to take picture of
        :param identify_exposure_time: exposure time of camera
        """
        ### acquire a single image
        self.camera.acquisition_stop()
        if identify_exposure_time is not None:
            initial_exposure_time = self.get_exposure_time()
            self.set_exposure_time(identify_exposure_time)
        if image_region is not None:
            initial_image_region = self.camera.image_region_get()
            self.camera.image_region_set(*image_region)
        self.camera.mode_shutter("Open")
        self.camera.mode_acquisition("Single Scan")
        self.camera.acquisition_start()
        self.wait_for_acquisition()

        ### return ANDOR to previous settings
        if image_region is not None:
            self.camera.image_region_set(*initial_image_region)
        if identify_exposure_time is not None:
            self.camera.setup_exposure_time(initial_exposure_time)

        self.camera.acquisition_stop()

    @rpc
    def continually_acquire_images(self, image_region=None, identify_exposure_time: TFloat=None) -> TNone:
        """
        Acquire many images from the camera
        todo: type annotation for arguments
        :param image_region: area of ion trap to take picture of
        :param identify_exposure_time: exposure time of camera
        """
        ### acquire many images
        self.camera.acquisition_stop()
        if identify_exposure_time is not None:
            self.set_exposure_time(identify_exposure_time)
        if image_region is not None:
            self.camera.image_region_set(*image_region)

        self.camera.mode_shutter("Open")
        self.camera.mode_acquisition("Run till abort")
        self.start_acquisition()
        self.camera.polling(True, 1.5)

    @rpc
    def stop_acquisition(self) -> TNone:
        """
        Stop camera from taking new images
        """
        self.camera.acquisition_stop()

    @rpc
    def start_acquisition(self) -> TNone:
        """
        Let Camera take new image(s)
        """
        self.camera.acquisition_start()

    @rpc
    def wait_for_acquisition(self) -> TNone:
        """
        Wait until image is acquired.
        Note: this function blocks until an image is acquired.
        """
        self.camera.acquisition_wait()


    def get_buffer_size(self) -> TInt64:
        """
        Get the maximum number of images that the circular buffer can store
        based on the current settings.
        Returns:
            int: the maximum number of images that can be stored.
        """
        return self.camera.buffer_size()

    @rpc
    def get_most_recent_image(self) -> TArray(TInt32, 1):
        """
        Retrieve most recent image camera took
        :return: TArray(TFloat): the most recent image the camera took
        """
        return self.camera.acquire_image_recent()

    @rpc
    def get_oldest_image(self,
                         notify_other_listeners: TBool = True) -> TArray(TInt32, 1):
        """
        Retrieve most oldest image camera took
        :return: TArray(TFloat): the most oldest image the camera took
        """
        return self.camera.acquire_image_oldest(notify_other_listeners)

    @rpc
    def get_oldest_images(self, num_images: TInt32):
        return self.camera.acquire_images_oldest(num_images)

    @rpc
    def get_all_acquired_images(self) -> TArray(TFloat, 1):
        """
        Retrieve all images in camera's data buffer
        :return: TArray(TFloat): the images in camera's data buffer
        """
        self.stop_acquisition()
        return self.camera.acquire_data()

    @rpc
    def set_exposure_time(self,
                          exposure_time: TFloat) -> TNone:
        """
        Set exposure time (in seconds) for image acquisition
        :param exposure_time: the exposure time (in seconds) to set.
        """

        self.camera.setup_exposure_time(exposure_time)

    @rpc
    def get_exposure_time(self) -> TFloat:
        """
        Get exposure time (in seconds) for image acquisition
        :return: exposure time in seconds
        """
        return self.camera.setup_exposure_time()

    @rpc
    def get_image_region(self) -> TTuple:
        """
        Get image region
        :return: tuple(horizontal bin, vertical bin, start x, stop x, start y, stop y) coordinates of image
        """
        return self.camera.image_region_get()


    @rpc
    def set_image_region(self, image_region: TTuple) -> TNone:
        self.camera.image_region_set(*image_region)

    @rpc
    def set_mode_shutter(self,
                         shutter_state):

        if shutter_state not in ['Open', 'Close', 'Auto']:
            raise ValueError('shutter state must be Open/Close/Auto')

        else:
            self.camera.mode_shutter(shutter_state)

    @rpc
    def set_mode_acquisition(self,
                         acquisition_state):

        if acquisition_state not in ['Single Scan', 'Accumulate', 'Kinetics', 'Fast Kinetics', 'Run till abort']:
            raise ValueError('acquistion state must be one of Single Scan/Accumulate/Kinetics/Fast Kinetics/Run till abort')

        else:
            self.camera.mode_acquisition(acquisition_state)

    @rpc
    def buffer_new_images(self) -> TTuple([TInt32, TInt32]):
        """
        Get the total number of new images in the buffer.
        Returns:
            (int, int): the indices of the first and last new image, respectively.
        """
        return self.camera.buffer_new_images()



