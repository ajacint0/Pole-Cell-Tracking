from skimage.measure import label, regionprops
import numpy as np
import tifffile as tfl
import csv
from glob import glob
import re
import os

os.chdir('/mnt/home/hnunley/ceph/tracked_embryos/2025-07-11_121341/gui_segs/')
for f in sorted(glob('*.tif')):
	num = re.findall('\d+', f)
	num0 = int(num[0])
	num1 = int(num[1])

	img = tfl.imread(f)
	nuclei_list = []
	nuclei = regionprops(img)
	for nucleus in nuclei:
		nuclei_list.append(f'{num0:03}_{nucleus.label:03}')
	print(nuclei_list)
	with open(f'/mnt/home/ajacinto/ceph/tracked_embryos/2025-07-11_121341/csvs/graph.csv', 'a', newline='') as f:
				writer = csv.writer(f)
				for nuc in nuclei_list:
				
					writer.writerow([nuc,'',''])
