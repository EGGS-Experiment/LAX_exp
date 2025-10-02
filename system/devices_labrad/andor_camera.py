from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad
from os import environ


class AndorCamera(LAXDevice):
    """
    High-level API functions for using the Andor Camera
    """
    name = "camera"

    def prepare_device(self):
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.camera = self.cxn.andor_server

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
            self.camera.set_image_region(*image_region)

        self.camera.mode_shutter("Open")
        self.camera.mode_acquisition("Run till abort")
        self.start_acquisition()
        self.camera.polling(True, 1.5)

    @rpc
    def stop_acquisition(self):
        """
        Stop camera from taking new images
        """
        self.camera.acquisition_stop()

    @rpc
    def start_acquisition(self):
        """
        Let Camera take new image(s)
        """
        self.camera.acquisition_start()

    @rpc
    def wait_for_acquisition(self):
        """
        Wait until image is acquired
        """
        self.camera.acquisition_wait()

    @rpc
    def get_most_recent_image(self) -> TArray(TInt32, 1):
        """
        Retrieve most recent image camera took
        :return: TArray(TFloat): the most recent image the camera took
        """
        return self.camera.acquire_image_recent()

    @rpc
    def get_all_acquired_images(self) -> TArray(TFloat, 1):
        """
        Retrieve all images in camera's data buffer
        :return: TArray(TFloat): the images in camera's data buffer
        """
        self.stop_acquisition()
        data = self.camera.acquire_data()
        return data

    @rpc
    def set_exposure_time(self, exposure_time: TFloat):
        """
        Set exposure time (in seconds) for image acquisition
        :param exposure_time: the exposure time (in seconds) to set.
        """
        self.stop_acquisition()
        self.camera.setup_exposure_time(exposure_time)
        self.start_acquisition()

    @rpc
    def get_exposure_time(self) -> TFloat:
        """
        Get exposure time (in seconds) for image acquisition
        :return: exposure time in seconds
        """
        self.stop_acquisition()
        return self.camera.setup_exposure_time()

    @rpc
    def set_image_region(self, image_region: TTuple) -> TNone:
        """
        Set image region
        :return: image_region (tuple): (horizontal bin, vertical bin, start x, stop x, start y, stop y) coordinates of image
        """
        self.stop_acquisition()
        self.camera.image_region_set(*image_region)
        self.start_acquisition()

    @rpc
    def get_image_region(self) -> TTuple:
        """
        Get image region
        :return: image_region (tuple): (horizontal bin, vertical bin, start x, stop x, start y, stop y) coordinates of image
        """
        self.stop_acquisition()
        return self.camera.image_region_get()
        self.start_acquisition()

