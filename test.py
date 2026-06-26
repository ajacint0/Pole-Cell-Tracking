from skimage.measure import label, regionprops
import numpy as np
import tifffile as tfl



img = tfl.imread('/mnt/home/ajacinto/ceph/tracked_embryos/2025-07-09_142109/cellpose/cropped_raw_ch_0_tp_42_cp_masks.tif')
props = regionprops(img)
for nucleus in props:
	if nucleus.label == 1047:
		print(nucleus.centroid)
#arr = np.load('/mnt/home/hnunley/tp_31_embryo_03_26_fixed.npy'

#print(arr)

'''
arr = np.load('/mnt/ceph/users/ajacinto/tracked_embryos/2025-03-26_144341/matches/new/matches_39_store_matches_new.npy')
print(arr)'''
