#DO THIS, SPLIT NUCLEUS 961 in TP 40

from skimage.measure import label, regionprops
import numpy as np
import tifffile as tfl
import networkx as nx
import matplotlib.pyplot as plt

def createCircularMask(shape, center, radius):

   grid = np.ogrid[tuple(slice(dim) for dim in shape)]
    
   # Create an array of coordinates with respect to the center
   coords = np.array([grid[i] - center[i] for i in range(len(shape))])

   # Calculate the distance of each point from the center of the sphere
   distances = np.sqrt(np.sum(coords ** 2, axis=0))
    
   # Create a mask based on the distance
   mask = distances <= radius

   return mask 

user = 'ajacinto'
movie = '2025-12-17_160948'
split_nucleus = False
combine_nuclei = False
add_nucleus = False
delete_nucleus = True
tp = 40
label = 34



img = tfl.imread(f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/cellpose/cropped_raw_ch_0_tp_{tp}_cp_masks.tif')
new_label1 = max(np.unique(img)) + 100
if split_nucleus == True:
	cut_z = False
	cut_x = True
	cutoffz = 218
	cutoffx = 222
	if cut_z == True:
		for z in range(0, img.shape[0]):
			if z < cutoffz:
				this_slice = img[z, :, :]
				this_slice[np.where(this_slice == label)] = new_label1
	elif cut_x == True:
		for x in range(0, img.shape[2]):
			if x < cutoffx:
				this_slice = img[:, :, x]
				this_slice[np.where(this_slice == label)] = new_label1

elif combine_nuclei == True:
	false_label = 35
	img[img == false_label] = label

elif add_nucleus == True:
	position = [213, 186, 240]
	mask = createCircularMask(img.shape, position, 10)
	img[mask] = new_label1
elif delete_nucleus == True:
	img[img == label] = 0


tfl.imwrite(f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/cellpose/cropped_raw_ch_0_tp_{tp}_cp_masks_new.tif' ,img)

