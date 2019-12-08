from __future__ import print_function
from skimage.segmentation import random_walker_segmentation
from skimage.segmentation import random_walker
import numpy as np
import nibabel as nib
import pydicom
from skimage import io
from skimage.external.tifffile import imshow
from skimage.transform import resize
from matplotlib import pyplot as plt
import time

def random_walker2d(image_path, show_image=True, save_path='random_walker2d.png'):
    """random walker 2d demo, support gray image and rgb 3 channel image

    # Args:
        image_path: str, image path
        show_image: bool, whether or not showing image
        save_path: bool, save the image pyplot show
    # Return:
        probility map, range in [0, 1.0]
    """
    img = io.imread(image_path)
    labels = np.zeros_like(img)
    # set walker seeds
    labels[0,:] = 1
    labels[-1,:] = 2
    now = time.time()
    conf_map = random_walker(img, labels, beta=130, mode='cg_mg', return_full_prob=True)
    print("Runtime", time.time() - now)
    if show_image:
        plt.figure('random walker 2d')
        plt.subplot(2, 1, 1)
        plt.imshow(img, cmap='gray')
        plt.subplot(2, 1, 2)
        plt.imshow(conf_map[0, :, :], cmap='gray')
        plt.savefig(save_path)

    return conf_map


def random_walker3d(dcm_path, target_shape=(200, 200, 400), show_image=True, 
                    save_nii=False, save_path='random_walker3d.nii'):
    """random walker 3d demo

    # Args:
        dcm_path: str, dicom file path
    # Return:
        probility map, range in [0, 1.0]
    """
    volume_dcm = pydicom.read_file(dcm_path)
    volume = volume_dcm.pixel_array
    original_shape = volume.shape
    print('orginal_volume shape is ', original_shape)
    volume = resize(volume, target_shape)
    print('resize volume to ', target_shape)
    labels = np.zeros_like(volume)
    # set random walker seeds
    labels[:, 0, :] = 2     
    labels[:, -1, :] = 1
    now = time.time()
    conf_map = random_walker(volume, labels, beta=90, mode='cg_mg', return_full_prob=True)
    conf_map = conf_map[1,:,:,:]
    print("运行时间", time.time() - now)
    if show_image:
        fig = plt.figure('random walker 3d')
        imshow(volume, figure=fig, subplot=(211), cmap='gray')
        imshow(conf_map, figure=fig, subplot=(212), cmap='gray')
        plt.show()
    if save_nii:
        conf_map = resize(conf_map, original_shape, preserve_range=True)
        conf_map = (conf_map * 255).astype('uint8')
        nii_volume = nib.Nifti1Image(conf_map, np.eye(4))
        nib.save(nii_volume, save_path)


if __name__ == '__main__':
    # random_walker2d('1.bmp')
    random_walker3d('./images/test.dcm')