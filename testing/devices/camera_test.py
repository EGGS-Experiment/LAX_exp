import numpy as np
from artiq.experiment import *


from LAX_exp.base import LAXExperiment
import cv2
from skimage import measure


class CameraTest(LAXExperiment, Experiment):
    """
    Experiment: Camera Test

    Test Camera
    """
    name = 'Camera Test'


    def build_experiment(self):

        # image region parameters: MAX (450,450) TO PREVENT LASER SCATTER OFF ELECTRODES FROM CONFUSING ANALYSIS
        self.setattr_argument('image_width_pixels', NumberValue(default=400,
                                                                     min=100, max=450, step=1,
                                                                     scale=1, ndecimals=0), group='Camera')

        self.setattr_argument('image_height_pixels', NumberValue(default=400,
                                                                   min=100, max=450, step=1,
                                                                   scale=1, ndecimals=0), group='Camera')

        self.setattr_device('camera')

    def prepare_experiment(self):
        self.image_region = (self.horizontal_image_region, self.vertical_image_region)

    @property
    def results_shape(self):
        return (2,2)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()
        self.camera.set_image_region(tuple(self.image_region))
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):

        self.camera.acquire_single_image()
        image = self.camera.get_most_recent_image()
        num_ions = self.show_ions(image)
        print(num_ions)

    # ANALYSIS
    def analyze_experiment(self):
        pass

    @rpc
    def show_ions(data) -> TFloat:


        plt.figure(1)
        plt.title("Original Image")
        plt.imshow(data)

        upper_percentile = np.percentile(data, 99)
        data[data < upper_percentile] = 0
        data[data >= upper_percentile] = 1

        kernel = np.ones((2, 2), np.uint8)
        data = cv2.erode(data, kernel, iterations=3)
        data = cv2.dilate(data, kernel, iterations=3)
        data[data > 0] = 1
        labels = measure.label(data)
        plt.figure(2)
        plt.imshow(data)
        plt.title("Manipulated Image")
        return len(np.unique(labels)) - 1



