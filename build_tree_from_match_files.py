import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import tifffile as tif
import csv


split_tp = 35
timevect = np.arange(40,113)

def createCircularMask(shape, center, radius):

   grid = np.ogrid[tuple(slice(dim) for dim in shape)]
    
   # Create an array of coordinates with respect to the center
   coords = np.array([grid[i] - center[i] for i in range(len(shape))])

   # Calculate the distance of each point from the center of the sphere
   distances = np.sqrt(np.sum(coords ** 2, axis=0))
    
   # Create a mask based on the distance
   mask = distances <= radius

   return mask 

movie = '2025-12-17_160948'
path_to_raw = f'/mnt/home/ajacinto/ceph/tracked_embryos/{movie}/test_cellpose/cropped_raw_ch_0_tp_'
path_to_nuclear_seg = f'/mnt/ceph/users/ajacinto/tracked_embryos/{movie}/gui_segs/tp_'

G = nx.Graph()
for i in timevect:


	npy_pairs = np.load(f'/mnt/home/ajacinto/ceph/tracked_embryos/{movie}/matches/matches_{i}_store_matches.npy')
	npy_pairs_scores = np.load(f'/mnt/home/ajacinto/ceph/tracked_embryos/{movie}/matches/matches_{i}_store_match_scores.npy')

	print(npy_pairs)

	for match_index in range(0, npy_pairs.shape[0]):
		if i == split_tp:
			assert npy_pairs[match_index, 0] <= 999
		elif i < split_tp:
			assert npy_pairs[match_index, 0] <= 999
			assert npy_pairs[match_index, 1] <= 999
		

		if npy_pairs_scores[match_index] > 0:
			if i == split_tp:
				G.add_edge(f'{i:03}_{npy_pairs[match_index, 0]:03}',f'{i+1:03}_{npy_pairs[match_index, 1]:06}')
			elif i > split_tp:
				G.add_edge(f'{i:03}_{npy_pairs[match_index, 0]:06}',f'{i+1:03}_{npy_pairs[match_index, 1]:06}')
			else:
				G.add_edge(f'{i:03}_{npy_pairs[match_index, 0]:03}',f'{i+1:03}_{npy_pairs[match_index, 1]:03}')
	
	node_count = 0
	for node in list(G.nodes):
		if int(node.split('_')[0]) == i:
			node_count = node_count + 1
	print(f'# of tp {i} nodes: {node_count}')

time_vals_vect = np.zeros_like(G.nodes)
labels_vals_vect = np.zeros_like(G.nodes)


for i in range(0, len(G.nodes)):
	nodename = list(G.nodes)[i]
	time_label = nodename.split('_')
	time_vals_vect[i] = int(time_label[0])
	labels_vals_vect[i] = int(time_label[1])

nuclear_intensity_vect = np.zeros_like(G.nodes)
nuclear_volume_vect = np.zeros_like(G.nodes)
nuclear_radius_val = 6
for timeval in timevect:
	V_raw = tif.imread(f'{path_to_raw}{timeval}.tif')
	try:
		file_path_to_nuclear_seg = f'{path_to_nuclear_seg}{timeval}_ch_1_seg.tif'
		Vnuclear_1 = tif.imread(file_path_to_nuclear_seg)
	except:
		try:
			file_path_to_nuclear_seg = f'{path_to_nuclear_seg}{timeval}_ch_1_seg_missing.tif'
			Vnuclear_1 = tif.imread(file_path_to_nuclear_seg)
		except:
			file_path_to_nuclear_seg = f'/mnt/home/ajacinto/ceph/tracked_embryos/{movie}/cellpose/cropped_raw_ch_0_tp_{timeval}_cp_masks.tif'
			Vnuclear_1 = tif.imread(file_path_to_nuclear_seg)

	unique_nuclear_labels = np.unique(Vnuclear_1)

	
	unique_labels = [int(i.split('_')[1]) for i in list(G.nodes) if int(i.split('_')[0]) == timeval]


	for i in unique_labels:
		if i == 0:
			continue
		nucleus = np.where(Vnuclear_1 == i)
		assert len(nucleus[0]) != 0
		try:
			if timeval > 35:
				findnodeval = list(G.nodes).index(f'{timeval:03}_{i:06}')
			else:
				findnodeval = list(G.nodes).index(f'{timeval:03}_{i:03}')
		except:
			continue

		mean_Z = int(np.mean(nucleus[0]))
		mean_Y = int(np.mean(nucleus[1]))
		mean_X = int(np.mean(nucleus[2]))
		
		volume = np.sum(Vnuclear_1 == i)

		nuclear_volume_vect[findnodeval] = volume

		mask = createCircularMask(Vnuclear_1.shape, [mean_Z, mean_Y, mean_X], nuclear_radius_val)
		
		nuclear_intensity_vect[findnodeval] = np.mean(V_raw[mask])
						
try:
	in_graph_file = []
	with open(f'/mnt/home/ajacinto/ceph/tracked_embryos/{movie}/csvs/volume_graph.csv', newline='') as file:
		csvfile = csv.reader(file)
		for lines in csvfile:
			in_graph_file.append(lines)

	print('huh')
	with open(f'/mnt/home/ajacinto/ceph/tracked_embryos/{movie}/csvs/volume_graph.csv', 'a') as csv_file:
		for i in range(0, len(list(G.nodes))):
			
			if [f'{list(G.nodes)[i]}',f'{nuclear_intensity_vect[i]}',f'{nuclear_volume_vect[i]}'] in in_graph_file:
			
				continue
			else:
				csv_file.write(list(G.nodes)[i])
				csv_file.write(',')
				csv_file.write(nuclear_intensity_vect[i])
				csv_file.write(',')
				csv_file.write(nuclear_volume_vect[i])
				csv_file.write('\n')
				
	csv_file.close()
except:
	pass
	
	with open(f'/mnt/home/ajacinto/ceph/tracked_embryos/{movie}/csvs/volume_graph.csv', 'w') as csv_file:
		for i in range(0, len(list(G.nodes))):
			csv_file.write(list(G.nodes)[i])
			csv_file.write(',')
			csv_file.write(nuclear_intensity_vect[i])
			csv_file.write(',')
			csv_file.write(nuclear_volume_vect[i])
			csv_file.write('\n')
	csv_file.close()


try:
	in_edge_file = []
	with open(f'/mnt/home/ajacinto/ceph/tracked_embryos/{movie}/csvs/volume_edges.csv', newline='') as file:
		csvfile = csv.reader(file)
		for lines in csvfile:
			print(f'adding {lines} to it')
			in_edge_file.append(lines)
	with open(f'/mnt/home/ajacinto/ceph/tracked_embryos/{movie}/csvs/volume_edges.csv', 'a') as csv_file:
		for i in range(0, len(list(G.edges))):
			
			if [f'{list(G.edges)[i][0]}',f'{list(G.edges)[i][1]}'] in in_edge_file:
			
				print(f'{list(G.edges)[i][0]},{list(G.edges)[i][1]} is in the file')
				continue
			else:
				print(f'{list(G.edges)[i][0]},{list(G.edges)[i][1]} is NOT in the file')
				csv_file.write(list(G.edges)[i][0])
				csv_file.write(',')
				csv_file.write(list(G.edges)[i][1])
				csv_file.write('\n')
except:
	print('doesnt exist')
	with open(f'/mnt/home/ajacinto/ceph/tracked_embryos/{movie}/csvs/volume_edges.csv', 'w') as csv_file:
		for i in range(0, len(list(G.edges))):
			csv_file.write(list(G.edges)[i][0])
			csv_file.write(',')
			csv_file.write(list(G.edges)[i][1])
			csv_file.write('\n')
			pass
	csv_file.close()


print(len(G.nodes))
