from artiq.experiment import *

from LAX_exp.base import LAXDevice

from os import environ
import labrad


class AndorCamera(LAXDevice):
    """
    High-level api functions for using the Andor Camera
    """

    name = "camera"

    def prepare_device(self):
        self.cxn = labrad.connect(environ['LABRADHOST'],
                                  port=7682, tls_mode='off',
                                  username='', password='lab')
        self.camera = self.cxn.andor_server

    @kernel(flags={"fast-math"})
    def initialize_device(self):
        pass

    @rpc
    def close_camera_shutter(self) -> TNone:
        """
        Close the Camera Shutter
        """
        self.camera.modeShutter('Close')

    @rpc
    def open_camera_shutter(self) -> TNone:
        """
        Open the Camera Shutter
        """
        self.camera.modeShutter('Open')

    @rpc
    def acquire_single_image(self, image_region, identify_exposure_time: TFloat)-> TNone:
        """
        Acquire a single image from the camera and then reset to previous setting

        Args:
            image_region: area of ion trap to take picture of

        Returns:
            image (np.array): Image acquired
        """
        ### acquire a single image
        self.camera.accquisition_stop()
        initial_exposure_time = self.camera.setup_exposure_time()
        self.camera.setup_exposure_time(identify_exposure_time)
        initial_image_region = self.camera.get_image_region()
        self.camera.set_image_region(*image_region)
        self.camera.mode_shutter("Open")
        self.camera.mode_accquisition("Single Scan")
        self.camera.accquisition_start()
        self.camera.accquisition_wait()

        ### return ANDOR to previous settings
        self.camera.setup_exposure_time(initial_exposure_time)
        self.camera.set_image_region(*initial_image_region)
        self.camera.mode_accquisition("Run till abort")
        self.camera.accquisition_start()
        self.camera.accquisition_wait()

    @rpc
    def continually_acquire_images(self, image_region, identify_exposure_time: TFloat) -> TNone:
        """
        Acquire a single image from the camera and then reset to previous setting

        Args:
            image_region: area of ion trap to take picture of

        Returns:
            image (np.array): Image acquired
        """
        ### acquire a single image
        self.camera.accquisition_stop()
        self.camera.setup_exposure_time(identify_exposure_time)
        self.camera.set_image_region(*image_region)
        self.camera.mode_shutter("Open")
        self.camera.mode_accquisition("Run till abort")
        self.camera.accquisition_start()

    @rpc
    def stop_acquisition_of_images(self) -> TNone:
        """
        Stop camera from taking new images
        """
        self.camera.accquisition_stop()

    @rpc
    def get_most_recent_image(self) -> TNone:
        """
        Retrieve most recent image camera took
        """
        return self.camera.acquireImageRecent()

    @rpc
    def get_all_acquired_images(self) -> TNone:
        """
        Retrieve all images in camera's data buffer
        """
        return self.camera.acquireData()

    @rpc
    def set_exposure_time(self, exposure_time: TFloat) -> TNone:
        """
        Set exposure time (in seconds) for image acquisition
        """
        self.camera.setup_exposure_time(exposure_time)

    @rpc
    def get_exposure_time(self) -> TFloat:
        """
        Get exposure time (in seconds) for image acquisition
        """
        return self.camera.setup_exposure_time()

    @rpc
    def set_image_region(self, image_region) -> TNone:
        """
        Set image region

        Args:
            image_region (tuple): (horizontal bin, vertical bin, start x, stop x, start y, stop y) coordinates of image
        """
        self.camera.setImageRegion(*image_region)

    @rpc
    def get_image_region(self):
        """
        Get image region

        Returns:
            image_region (tuple): (horizontal bin, vertical bin, start x, stop x, start y, stop y) coordinates of image
        """
        return self.camera.getImageRegion()



