import numpy as np
from artiq.experiment import *

from LAX_exp.base import LAXExperiment
from skimage import measure
import skimage
from artiq.language.units import *

from matplotlib import pyplot as plt

class CameraTest(LAXExperiment, Experiment):
    """
    Experiment: Camera Test

    Test Camera
    """
    name = 'Camera Test'
    IMAGE_HEIGHT = 512
    IMAGE_WIDTH = 512


    def build_experiment(self):

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

        self.setattr_device('camera')

    def prepare_experiment(self):

        border_x = (self.IMAGE_WIDTH - self.image_width_pixels) / 2
        border_y = (self.IMAGE_HEIGHT - self.image_height_pixels) / 2

        start_x = int(border_x)
        end_x = int(self.IMAGE_WIDTH - border_x) - 1

        start_y = int(border_y)
        end_y = int(self.IMAGE_WIDTH - border_y) - 1

        print(ms)

        self.image_region = (self.horizontal_binning, self.vertical_binning,
                             start_x, end_x, start_y, end_y)

    @property
    def results_shape(self):
        return (2,2)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):


        self.core.break_realtime()
        self.camera.set_image_region(self.image_region)
        self.core.break_realtime()


    @kernel(flags={"fast-math"})
    def run_main(self):

        self.core.break_realtime()
        self.camera.continually_acquire_images()
        delay(5*s)
        self.core.break_realtime()
        data = self.camera.get_all_acquired_images()
        self.reshape_image(data)



        # self.camera.acquire_single_image()
        # image = self.camera.get_most_recent_image()
        # num_ions = self.show_ions(image)
        # print(num_ions)




    # ANALYSIS
    def analyze_experiment(self):
        pass

    @rpc
    def reshape_image(self, image):
        print(np.shape(np.reshape(image, (-1, self.image_width_pixels, self.image_width_pixels))))

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
        labels = measure.label(data)
        plt.figure(2)
        plt.imshow(data)
        plt.title("Manipulated Image")
        plt.savefig("Z:\motion\Pictures\manipulated.png")
        return len(np.unique(labels)) - 1



