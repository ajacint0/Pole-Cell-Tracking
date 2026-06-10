from cellpose import models
import tifffile as tfl
import numpy as np
import os
from glob import glob
import re

split_tp = 35
os.chdir('/mnt/home/ajacinto/ceph/tracked_embryos/2025-12-17_160948/new_test_cellpose/')
for f in sorted(glob('*.tif')):
	num = re.findall('\d+', f)
	num0 = int(num[0])
	num1 = int(num[1])
	if num1 < split_tp:
		continue

	img = tfl.imread(f)
	model = models.Cellpose(gpu=True, model_type='cyto2')
	masks, flows, styles, diams = model.eval(img, tile_overlap = .5, do_3D=True,anisotropy = 1.0, tile=True, min_size=128, batch_size=8)
	tfl.imwrite(f'/mnt/home/ajacinto/ceph/tracked_embryos/2025-12-17_160948/new_cellpose/cropped_raw_ch_{num0}_tp_{num1}_cp_masks.tif', masks.astype(np.uint16))
