from __future__ import print_function, unicode_literals, absolute_import, division
import numpy as np
import matplotlib
matplotlib.rcParams["image.interpolation"] = 'none'
import matplotlib.pyplot as plt
from glob import glob
import tifffile as tfl
from skimage.morphology import remove_small_objects
from skimage.measure import label, regionprops
from skimage.transform import rescale, warp
from skimage.segmentation import find_boundaries
from scipy.ndimage import affine_transform
from point_cloud_registration import ICP, PlaneICP, NDT, VPlaneICP
import copy
import random
import csv

#matplotlib.use('TkAgg')


def register_volumes(movie, source_tp, source, target, matrix):
	if np.array_equal(matrix, np.eye(4)):

		source_segmentation = source
		source_nuclei = remove_small_objects(source_segmentation)
		

		try:
			source_outline = find_boundaries(source_segmentation,mode='inner')
			#source_outline = source_segmentation
			source_input = np.argwhere(source_outline > 0)
			source_indices = np.random.choice(source_input.shape[0], int(source_input.shape[0] * .05), replace = False)
			source_input = source_input[source_indices]


			target_segmentation = target
			target_nuclei = remove_small_objects(target_segmentation)

			target_outline = find_boundaries(target_segmentation,mode='inner')
			#target_outline = target_segmentation

			target_input = np.argwhere(target_outline > 0)
			target_indices = np.random.choice(target_input.shape[0], int(target_input.shape[0] * .05), replace = False)
			target_input = target_input[target_indices]


			icp = VPlaneICP(voxel_size = 10, max_iter=60, max_dist = 20, tol=1e-3)
			icp.set_target(target_input)
			T_new = icp.align(source_input, init_T=np.eye(4))
		except:
			
			source_outline = source_segmentation
			source_input = np.argwhere(source_outline > 0)
			source_indices = np.random.choice(source_input.shape[0], int(source_input.shape[0] * .05), replace = False)
			source_input = source_input[source_indices]


			target_segmentation = target
			target_nuclei = remove_small_objects(target_segmentation)

			
			target_outline = target_segmentation

			target_input = np.argwhere(target_outline > 0)
			target_indices = np.random.choice(target_input.shape[0], int(target_input.shape[0] * .05), replace = False)
			target_input = target_input[target_indices]


			icp = VPlaneICP(voxel_size = 10, max_iter=60, max_dist = 20, tol=1e-3)
			icp.set_target(target_input)
			T_new = icp.align(source_input, init_T=np.eye(4))
		T_new = (np.array(T_new, dtype=np.float64))
		np.save(f'/mnt/home/ajacinto/ceph/tracked_embryos/{movie}/transformations/tp_{source_tp}_transformation.npy', T_new)
	else:
		T_new = matrix
	#ones = np.ones((end_source.shape[0], 1))
	#points_h = np.hstack([end_source, ones])

	#transformed_points_h = points_h @ T_new.T   
	#transformed_points = transformed_points_h[:, :3]
	source_segmentation = source
	rotation = T_new[:3, :3]
	translation = T_new[:3, 3]
	
	inv_rotation = np.linalg.inv(rotation)
	inv_translation = -inv_rotation @ translation

	transformed_image = affine_transform(source_segmentation, matrix = inv_rotation, offset = inv_translation, order=0)
	
	
	
	return transformed_image

