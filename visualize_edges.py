import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

import networkx as nx
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import tifffile as tif
import csv
from skimage.measure import regionprops, label
from skimage.morphology import remove_small_objects
from trim_graph import trim_graph
from register_volumes import register_volumes


register = False
split_tp = 35
G = nx.Graph()
cutoff = 0
limit_graph = True

hayden_arr = True
highlight = []
text = True




path_to_nuclear_seg = '/mnt/ceph/users/ajacinto/tracked_embryos/2025-12-17_160948/gui_segs/tp_'
path_to_membrane_seg = '/mnt/home/ajacinto/ceph/tracked_embryos/2025-12-17_160948/cellpose/cropped_raw_ch_0_tp_'



# Loads nuclei that will go in lineage tree
if hayden_arr == True:
	maybe_arr = np.load('/mnt/home/ajacinto/ceph/tracked_embryos/2025-12-17_160948/important_nuclei/tp_35_maybe.npy')
	for_sure_arr = np.load('/mnt/home/ajacinto/ceph/tracked_embryos/2025-12-17_160948/important_nuclei/tp_35_for_sure.npy')
	
	arr = []
	for i in maybe_arr:
		arr.append(i)
	for i in for_sure_arr:
		arr.append(i)
	print(arr)

#Chooses which timepoints are shown on screen, can be far apart like tp_early=4, tp_late=30
tp_early = 54
tp_late = 55


#split_tp is when I switch from manual segmentation to cellpose segmentation, with cellpose segmentations have thousands of labels, I need to change the naming convention(I should change this to all be 6 characters long, made more sense when manually writing onto csv)
if tp_early <= split_tp:
	early_buffer = 3
else:
	early_buffer = 6

if tp_late <= split_tp:
	late_buffer = 3
else:
	late_buffer = 6

start_edge = []
end_edge = []


#Adds nodes and edges to graph
with open('/mnt/home/ajacinto/ceph/tracked_embryos/2025-12-17_160948/csvs/volume_edges.csv', newline='') as file:
	csvfile = csv.reader(file)
	data = list(csvfile)
	for lines in data:
		G.add_edge(lines[0],lines[1])
		

only_these_nuclei = []

#Adds nodes with no edges
with open('/mnt/home/ajacinto/ceph/tracked_embryos/2025-12-17_160948/csvs/volume_graph.csv', newline='') as file:
		csvfile = csv.reader(file)
		for lines in csvfile:
			G.add_node(lines[0])


#If false, shows every segmentation that is large enough in size
if limit_graph == False:

	with open('/mnt/home/ajacinto/ceph/tracked_embryos/2025-12-17_160948/csvs/volume_graph.csv', newline='') as file:
		csvfile = csv.reader(file)
		for lines in csvfile:
			if f'{tp_early:03}' in (lines[0])[:3]:
				only_these_nuclei.append(lines[0])
				start_edge.append(lines[0])
			if f'{tp_late:03}' in (lines[0])[:3]:
				only_these_nuclei.append(lines[0])
				end_edge.append(lines[0])

#Else, if using the npy file containing important nuclei, save them for highlighting with bright colors, if not using npy file, choose cutoff value, only highlights nuclei that have more connections in the lineage than the cutoff value
else:

	if hayden_arr == True:
		for node in arr:
			#trim_graph tells the program which tracks to highlight based on the minimum number of nodes in each track
			node_list, delete_list = trim_graph(0, node, data, [node], [], cutoff)
			for vertex in node_list:
				try:

					highlight.append(vertex) 
				except:
					continue

	else:

		for i in range(1, 83):
			
			for node in list(G):
				if int(node.split('_')[0]) != i:
					continue
				
				skip = False
				for edge in data:
					if node == edge[1]:
						skip = True
						break
				if skip == True:
					continue
				else:

					node_list, delete_list = trim_graph(0, node, data, [node], [], cutoff)

				
					if False in delete_list:
						for vertex in node_list:
							try:
									
								highlight.append(vertex)
							except:
								continue

	
	for node in list(G.nodes):
		if f'{tp_early:03}' in node[:3]:
			only_these_nuclei.append(node)
			start_edge.append(node)
		if f'{tp_late:03}' in node[:3]:
			if node == '035_001286':
				print('ahah')
			only_these_nuclei.append(node)
			end_edge.append(node)

#Reads in the sementations
print(highlight)
try:
	file_path_to_nuclear_seg_early = f'{path_to_nuclear_seg}{tp_early}_ch_1_seg.tif'
	nuclear_seg_early = tif.imread(file_path_to_nuclear_seg_early)
except FileNotFoundError:
	try:
		file_path_to_nuclear_seg_early = f'{path_to_nuclear_seg}{tp_early}_ch_1_seg_missing.tif'
		nuclear_seg_early = tif.imread(file_path_to_nuclear_seg_early)
	except FileNotFoundError:
		file_path_to_nuclear_seg_early = f'{path_to_membrane_seg}{tp_early}_cp_masks.tif'
		nuclear_seg_early = tif.imread(file_path_to_nuclear_seg_early)
		nuclear_seg_early = remove_small_objects(nuclear_seg_early)
try:
	file_path_to_nuclear_seg_late = f'{path_to_nuclear_seg}{tp_late}_ch_1_seg.tif'
	nuclear_seg_late = tif.imread(file_path_to_nuclear_seg_late)
except FileNotFoundError:
	try:
		file_path_to_nuclear_seg_late = f'{path_to_nuclear_seg}{tp_late}_ch_1_seg_missing.tif'
		nuclear_seg_late = tif.imread(file_path_to_nuclear_seg_late)
	except FileNotFoundError:
		file_path_to_nuclear_seg_late = f'{path_to_membrane_seg}{tp_late}_cp_masks.tif'
		nuclear_seg_late = tif.imread(file_path_to_nuclear_seg_late)
		nuclear_seg_late = remove_small_objects(nuclear_seg_late)

#If true, the window that show's the graph will have the centroids registered (they were already registered when doing IOU calculations to make the initial edges)
if register == True:
	matrix = np.load(f'/mnt/home/ajacinto/ceph/tracked_embryos/2025-12-17_160948/transformations/tp_{tp_early}_transformation.npy')
	nuclear_seg_early = register_volumes('2025-12-17_160948', tp_early, nuclear_seg_early, nuclear_seg_late, matrix)

props_early = regionprops(nuclear_seg_early)
props_late = regionprops(nuclear_seg_late)


x_early = []
y_early = []
z_early = []

x_late = []
y_late = []
z_late = []

color_early = []
color_late = []

name_early = []
name_late = []

fig = plt.figure()
ax = plt.axes(projection='3d')

#gets nucleus centroid coords, assigns color, assigns text to nodes in graph for early timepoint

for nucleus in props_early:
	print(f'nucleus {nucleus.label}: {nucleus.centroid}')
	#print(nucleus.label)
	if f'{tp_early:03}_{nucleus.label:0{early_buffer}}' not in only_these_nuclei and early_buffer == 3:
		ax.scatter(nucleus.centroid[2],nucleus.centroid[1],nucleus.centroid[0], s=10, c='brown', marker='o')
		if text == True:
			ax.text(nucleus.centroid[2], nucleus.centroid[1], nucleus.centroid[0], f'{tp_early:03}_{nucleus.label:0{early_buffer}}')
		continue
	if f'{tp_early:03}_{nucleus.label:0{early_buffer}}' in only_these_nuclei:
		if f'{tp_early:03}_{nucleus.label:0{early_buffer}}' in highlight:
			x_early.append(nucleus.centroid[2])
			y_early.append(nucleus.centroid[1])
			z_early.append(nucleus.centroid[0])
			color_early.append('yellow')
			name_early.append(f'{tp_early:03}_{nucleus.label:0{early_buffer}}')
			if text == True:
				ax.text(nucleus.centroid[2], nucleus.centroid[1], nucleus.centroid[0], f'{tp_early:03}_{nucleus.label:0{early_buffer}}')
		else:
			x_early.append(nucleus.centroid[2])
			y_early.append(nucleus.centroid[1])
			z_early.append(nucleus.centroid[0])
			color_early.append('red')
			name_early.append(f'{tp_early:03}_{nucleus.label:0{early_buffer}}')
			if text == True:
				ax.text(nucleus.centroid[2], nucleus.centroid[1], nucleus.centroid[0], f'{tp_early:03}_{nucleus.label:0{early_buffer}}')

ax.scatter(x_early,y_early,z_early, s=10, c=color_early, marker='o')

#gets nucleus centroid coords, assigns color, assigns text to nodes in graph for late timepoint

for nucleus in props_late:
	
	if f'{tp_late:03}_{nucleus.label:0{late_buffer}}' not in only_these_nuclei and late_buffer == 3:

		ax.scatter(nucleus.centroid[2],nucleus.centroid[1],nucleus.centroid[0], s=10, c='green', marker='o')
		if text == True:
			ax.text(nucleus.centroid[2], nucleus.centroid[1], nucleus.centroid[0], f'{tp_late:03}_{nucleus.label:0{late_buffer}}')
		continue
	if f'{tp_late:03}_{nucleus.label:0{late_buffer}}' in only_these_nuclei:
		if f'{tp_late:03}_{nucleus.label:0{late_buffer}}' in highlight:
			x_late.append(nucleus.centroid[2])
			y_late.append(nucleus.centroid[1])
			z_late.append(nucleus.centroid[0])
			color_late.append('lime')
			name_late.append(f'{tp_late:03}_{nucleus.label:0{late_buffer}}')
			if text == True:
				ax.text(nucleus.centroid[2], nucleus.centroid[1], nucleus.centroid[0], f'{tp_late:03}_{nucleus.label:0{late_buffer}}')
		else:
			x_late.append(nucleus.centroid[2])
			y_late.append(nucleus.centroid[1])
			z_late.append(nucleus.centroid[0])
			color_late.append('blue')
			name_late.append(f'{tp_late:03}_{nucleus.label:0{late_buffer}}')
			if text == True:
				ax.text(nucleus.centroid[2], nucleus.centroid[1], nucleus.centroid[0], f'{tp_late:03}_{nucleus.label:0{late_buffer}}')

ax.scatter(x_late,y_late,z_late, s=10, c=color_late, marker='o')


#finds a path between the early and late nuclei even if they are multiple timepoints apart

edges_to_show = []

for i in range(0, len(start_edge)):
	for j in range(0, len(end_edge)):

		try:
			distance, path = nx.single_source_dijkstra(G, start_edge[i], end_edge[j])
			add_edge = True
			for k in range(0, len(path) - 1):
				if int(path[k].split('_')[0]) > int(path[k + 1].split('_')[0]):
					add_edge = False

			if add_edge == True:
				edges_to_show.append([path[0], path[-1]])
		except:
			pass

#plots the text for the graph

for edge in edges_to_show:
	ind_early = name_early.index(edge[0])
	
	ind_late = name_late.index(edge[1])

	ax.plot([x_early[ind_early], x_late[ind_late]],[y_early[ind_early], y_late[ind_late]],[z_early[ind_early], z_late[ind_late]], color = 'white')
	
ax.set_facecolor('purple')
fig.patch.set_facecolor('black')


ax.set_xticks(range(0,(nuclear_seg_late.shape[2]), 100))
ax.set_yticks(range(0,(nuclear_seg_late.shape[1]), 30))
ax.set_zticks(range(0,(nuclear_seg_late.shape[0]), 100))
plt.show()











	
	








