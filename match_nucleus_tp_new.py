
from __future__ import print_function, unicode_literals, absolute_import, division
import numpy as np
import matplotlib
matplotlib.rcParams["image.interpolation"] = 'none'
import matplotlib.pyplot as plt
from glob import glob
#from tqdm import tqdm
from tifffile import imread
#from csbdeep.utils import Path, download_and_extract_zip_file
#from stardist import relabel_image_stardist3D, Rays_GoldenSpiral, calculate_extents
#from stardist import fill_label_holes, random_label_cmap
#from stardist.matching import matching_dataset, matching
from skimage.morphology import remove_small_objects
from scipy.optimize import linear_sum_assignment
from register_volumes import register_volumes
import os
#source /mnt/home/ajacinto/venvs/my_env/bin/activate

movie = '2025-12-17_160948'
user = 'ajacinto'
def iou_matching(img1, img2, timeval):
	iou_matrix = np.zeros((len(np.unique(img1)),len(np.unique(img2))))
	counter = 0
	for label1 in np.unique(img1):
		
		if label1 == 0:
			counter = counter + 1
			continue
		voxel_values = img2[np.where(img1 == label1)]
		unique_voxel_values = np.unique(voxel_values)
		for label2 in unique_voxel_values:
			print(label2)
			if label2 > 0:
				this_label2 = np.where(np.unique(img2) == label2)
				intersect_val = len(np.where(voxel_values == label2)[0])
				np_or = np.logical_or(img1 == label1, img2 == label2)
				np_or = np.where(np_or)
				union_value = len(np_or[0])
				iou_matrix[counter, this_label2] = intersect_val / union_value

		counter = counter + 1
	print(iou_matrix)
	M = linear_sum_assignment(iou_matrix, True)
	M_again = (np.zeros(len(M[0])),np.zeros(len(M[1])))
	print(M)
	matches = []
	match_scores = []
	for i in range(1, len(M[0])):
		M1 = M[0][i]
		M2 = M[1][i]
		M_again[0][i] = np.unique(img1)[M1]
		M_again[1][i] = np.unique(img2)[M2]
		print('here')
		print(M_again[0][i])
		print(M_again[1][i])
		matches.append([int(M_again[0][i]), int(M_again[1][i])])
		match_scores.append(iou_matrix[M[0][i]][M[1][i]])
	np.save(path_for_saving_matches + str(timeval) + '_store_matches_new.npy', matches)
	np.save(path_for_saving_matches + str(timeval) + '_store_match_scores_new.npy', match_scores)



dirs = ['matches', 'csvs', 'transformations']
for dir_ in dirs:
	dir_path = f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/{dir_}/'
	if not os.path.isdir(dir_path):
		os.makedirs(f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/{dir_}/', exist_ok=True)

split_tp = 35



timevect = np.arange(75,114) #20,27

path_to_nuclear_segmentations = f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/gui_segs/tp_'
path_to_membrane_segmentations = f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/cellpose/cropped_raw_ch_0_tp_'



path_for_saving_matches = f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/matches/matches_'

suffix_for_path = '_ch_1_seg.tif'
alt_suffix_for_path = '_ch_1_seg_missing.tif'

for timeval in timevect:
	if timeval + 1 <= split_tp:
		try:
			X1 = imread(path_to_nuclear_segmentations + str(timeval + 1) + suffix_for_path)
		except:
			X1 = imread(path_to_nuclear_segmentations + str(timeval + 1) + alt_suffix_for_path)
	else:
		X1 = imread(path_to_membrane_segmentations + str(timeval + 1) + '_cp_masks.tif')
		X1 = remove_small_objects(X1, 5000)

	if timeval <= split_tp:
		try:
			Xi = imread(path_to_nuclear_segmentations + str(timeval) + suffix_for_path)
		except:
			Xi = imread(path_to_nuclear_segmentations + str(timeval) + alt_suffix_for_path)
	else:
		Xi = imread(path_to_membrane_segmentations + str(timeval) + '_cp_masks.tif')
		Xi = remove_small_objects(Xi, 5000)

	print(f'tp: {timeval}')
	Xi_transformed = register_volumes(movie, timeval, Xi, X1, np.eye(4))
	iou_matching(Xi_transformed, X1, timeval)
