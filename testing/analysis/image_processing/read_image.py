import numpy as np
import cv2
import matplotlib.pyplot as plt
import skimage.morphology
from scipy.ndimage import gaussian_filter
from skimage import measure
import pandas as pd

"""
Inspiration from: https://pyimagesearch.com/2016/10/31/detecting-multiple-bright-spots-in-an-image-with-python-and-opencv/
"""

file_path = r'C:\Users\joshr\OneDrive\Documents\UCLA\Research\Example_Data\single_ion_7_20_Wed_Jul_20_2022_36.fits'
file_path = r'C:\Users\joshr\OneDrive\Documents\UCLA\Research\Example_Data\2024_05_02_image_signal(1).csv'
file_path = r'C:\Users\joshr\OneDrive\Documents\UCLA\Research\Example_Data\2024_05_02_image_signal_flippedish.csv'


def show_ions(data):

    plt.figure(1)
    plt.title("UnCropped Image")
    plt.imshow(data)

    if np.shape(data)[0] > 450:
        data = data[31:481, :]

    if np.shape(data)[1] > 450:
        data = data[:, 31:481]

    plt.figure(2)
    plt.title("Cropped Image")
    plt.imshow(data)


    upper_percentile = np.percentile(data, 99)
    data[data < upper_percentile] = 0
    data[data >= upper_percentile] = 1


    kernel = np.ones((2, 2), np.uint8)
    for i in range(3):
        data = np.uint16(skimage.morphology.binary_erosion(data, kernel))

    for i in range(3):
        data = np.uint16(skimage.morphology.binary_dilation(data, kernel))
    # data = cv2.dilate(data, kernel, iterations=3)
    data[data > 0] = 1
    data[data > 0] = 1
    labels = measure.label(data)
    plt.figure(3)
    plt.imshow(data)
    plt.title("Manipulated Image")
    return len(np.unique(labels))-1


from astropy.io import fits

# with fits.open(file_path) as hdul:  # open a FITS file
#
#     data = hdul[0].data  # assume the first extension is an imag
#     # data = np.array(data[0, 150:350, 150:350])
#     data = np.array(data[0,:,:])
#     data_org = data
#
# # thr
# data[data < np.percentile(data, 99)] = 0
# data[data >= np.percentile(data, 99)] = 1
# kernel = np.ones((2, 2), np.uint8)
# data = cv2.erode(data, kernel, iterations=2)
# data = cv2.dilate(data, kernel, iterations=3)
# data[data > 0] = 1
# labels = measure.label(data)
# print(len(np.unique(labels)) - 1)
#
# plt.figure(1)
# plt.imshow(data, cmap='gray')
# plt.figure(2)
# plt.imshow(data_org, cmap='gray')
# plt.show()

df = pd.read_csv(file_path, header=None)

image_vals = np.array(np.reshape(df, len(df)))
image = np.reshape(image_vals, (512, 512))
print(show_ions(image))
plt.show()
