#encoding: utf-8
from __future__ import print_function
from skimage import io
from confidence_map import confidence_map3d, confidence_map2d
import numpy as np
import pydicom
from skimage.external.tifffile import imshow
import time
from matplotlib import pyplot as plt
from skimage.exposure import rescale_intensity
from skimage.transform import resize


def conf3d_demo(path):
    volume_dcm = pydicom.read_file(path)
    volume = volume_dcm.pixel_array
    slice_spacing = volume_dcm.SpacingBetweenSlices # 0.520
    height_spacing = volume_dcm.PixelSpacing[0]     # 0.082
    width_spacing = volume_dcm.PixelSpacing[1]      # 0.2
    volume = np.transpose(volume, [1, 2, 0])        # depth last
    volume = volume[:,:,:20]
    v_min, v_max = (0.007, 0.81)
    volume = rescale_intensity(volume, in_range=(v_min * 255, v_max * 255))
    print(volume.shape)
    now = time.time()
    conf_map = confidence_map3d(volume, alpha=1.5, beta=90, gamma=0.03, solver_mode='gpu')
    print(conf_map.shape)
    print("Runtime ", time.time() - now)
    conf_map = np.transpose(conf_map, [2, 0, 1])
    volume = np.transpose(volume, [2, 0, 1])
    fig = plt.figure()
    imshow(conf_map, figure=fig, subplot=(211), cmap='gray')
    imshow(volume, figure=fig, subplot=(212), cmap='gray')
    plt.show()


def conf2d_demo(image_path, v_min=0.007, v_max=0.81, threshold=0.5):
    """confidence map 2d demo
    
    # Args
        image_path: str, image path
        v_min, v_max: scale image intensity to range of [v_min, v_max]

    # Retrun
        confience map, actully it is a probility map with range [0, 1.0]
    """
    img = io.imread(image_path)
    img = rescale_intensity(img, in_range=(v_min*255, v_max*255))
    height_spacing = 0.082
    width_spacing = 0.2
    print('image shape is ', img.shape)
    now = time.time()
    spacing = [1.0, width_spacing / height_spacing]
    conf_map = confidence_map2d(img, alpha=1.5, beta=90, gamma=0.03, spacing=None, solver_mode='bf')
    conf_map = np.clip(conf_map, 0, 1)
    print("Runtime ", time.time() - now)
    plt.subplot(221)
    plt.imshow(img)
    plt.axis('off')
    plt.subplot(223)
    plt.imshow(conf_map)
    plt.axis('off')
    plt.subplot(222)
    plt.imshow((conf_map > threshold).astype('uint8')*img)
    plt.axis('off')
    plt.subplot(224)
    plt.imshow((conf_map > threshold).astype('uint8'))
    plt.axis('off')
    plt.savefig('confidence2d.png')

    return conf_map


if __name__ == "__main__":
    # conf3d_demo('./test.dcm')
    conf2d_demo('./images/test.bmp')