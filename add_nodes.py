from skimage.measure import label, regionprops
import numpy as np
import tifffile as tfl
import csv

nuclear_seg = '/mnt/home/hnunley/ceph/tracked_embryos/2025-07-11_121341/gui_segs/tp_11_ch_1_seg.tif'
tp = 11
img = tfl.imread(nuclear_seg)
nuclei_list = []
nuclei = regionprops(img)
for nucleus in nuclei:
	nuclei_list.append(f'{tp:03}_{nucleus.label:03}')
print(nuclei_list)
with open(f'/mnt/home/ajacinto/ceph/tracked_embryos/2025-07-11_121341/csvs/graph.csv', 'a', newline='') as f:
			writer = csv.writer(f)
			for nuc in nuclei_list:
			
				writer.writerow([nuc,'',''])

