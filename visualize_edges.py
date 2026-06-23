import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backend_bases import MouseButton
import matplotlib.pyplot as plt
import mplcursors
import networkx as nx
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import tifffile as tif
import csv
from skimage.measure import regionprops, label
from skimage.morphology import remove_small_objects
from trim_graph import trim_graph
from register_volumes import register_volumes
from mpl_toolkits.mplot3d import proj3d
import sys

user = 'ajacinto'
movie = '2025-07-09_142109'

register = False
split_tp = 0
G = nx.Graph()
cutoff = 0
limit_graph = True

important_arr = True
highlight = []
text = True




path_to_nuclear_seg = f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/gui_segs/tp_'
path_to_membrane_seg = f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/cellpose/cropped_raw_ch_0_tp_'



# Loads nuclei that will go in lineage tree
if important_arr == True:
	maybe_arr = np.load(f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/important_nuclei/tp_6_maybe.npy')
	for_sure_arr = np.load(f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/important_nuclei/tp_6_for_sure.npy')
	
	arr = []
	for i in maybe_arr:
		arr.append(i)
	for i in for_sure_arr:
		arr.append(i)
	#print(arr)

#Chooses which timepoints are shown on screen, can be far apart like tp_early=4, tp_late=30
tp_early = 13
tp_late = 14


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
with open(f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/csvs/test/edges.csv', newline='') as file:
	csvfile = csv.reader(file)
	data = list(csvfile)
	for lines in data:
		G.add_edge(lines[0],lines[1])
		

only_these_nuclei = []

#Adds nodes with no edges
with open(f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/csvs/graph.csv', newline='') as file:
		csvfile = csv.reader(file)
		for lines in csvfile:
			G.add_node(lines[0])


#If false, shows every segmentation that is large enough in size
if limit_graph == False:

	with open(f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/csvs/graph.csv', newline='') as file:
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

	if important_arr == True:
		for node in arr:
			#trim_graph tells the program which tracks to highlight based on the minimum number of nodes in each track
			node_list, delete_list = trim_graph(0, node, data, [node], [], cutoff)
			for vertex in node_list:
				try:

					highlight.append(vertex) 
				except:
					continue

	else:

		for i in range(1, 56):
			
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
	matrix = np.load(f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/transformations/tp_{tp_early}_transformation.npy')
	nuclear_seg_early = register_volumes(movie, tp_early, nuclear_seg_early, nuclear_seg_late, matrix)

props_early = regionprops(nuclear_seg_early)
props_late = regionprops(nuclear_seg_late)


x = np.array([])
y = np.array([])
z = np.array([])

x_late = []
y_late = []
z_late = []

color = np.array([])
color_late = []

name = np.array([])
name_late = []

early_txt = []

fig = plt.figure()
ax = plt.axes(projection='3d')

#gets nucleus centroid coords, assigns color, assigns text to nodes in graph for early timepoint

for nucleus in props_early:
	#print(f'nucleus {nucleus.label}: {nucleus.centroid}')
	#print(nucleus.label)
	if f'{tp_early:03}_{nucleus.label:0{early_buffer}}' not in only_these_nuclei and early_buffer == 3:
		ax.scatter(nucleus.centroid[2],nucleus.centroid[1],nucleus.centroid[0], s=10, c='brown', marker='o')
		if text == True:
			ax.text(nucleus.centroid[2], nucleus.centroid[1], nucleus.centroid[0], f'{tp_early:03}_{nucleus.label:0{early_buffer}}')
		continue
	if f'{tp_early:03}_{nucleus.label:0{early_buffer}}' in only_these_nuclei:
		if f'{tp_early:03}_{nucleus.label:0{early_buffer}}' in highlight:
			#print('hi')
			#print(nucleus.centroid[2])
			x=np.append(x, nucleus.centroid[2])
			y=np.append(y,nucleus.centroid[1])
			z=np.append(z,nucleus.centroid[0])
			color=np.append(color,'yellow')
			name=np.append(name,f'{tp_early:03}_{nucleus.label:0{early_buffer}}')
			
			txt = ax.text(nucleus.centroid[2], nucleus.centroid[1], nucleus.centroid[0], f'{tp_early:03}_{nucleus.label:0{early_buffer}}')
			txt.set_visible(False)
			early_txt.append(txt)
				
		else:
			x= np.append(x,nucleus.centroid[2])
			y=np.append(y,nucleus.centroid[1])
			z=np.append(z,nucleus.centroid[0])
			color=np.append(color,'red')
			name=np.append(name,f'{tp_early:03}_{nucleus.label:0{early_buffer}}')
			
			txt = ax.text(nucleus.centroid[2], nucleus.centroid[1], nucleus.centroid[0], f'{tp_early:03}_{nucleus.label:0{early_buffer}}')
			txt.set_visible(False)
			early_txt.append(txt)

#sc_early = ax.scatter(x_early,y_early,z_early, s=10, c=color_early, marker='o')
#print(x)

#gets nucleus centroid coords, assigns color, assigns text to nodes in graph for late timepoint

for nucleus in props_late:
	
	if f'{tp_late:03}_{nucleus.label:0{late_buffer}}' not in only_these_nuclei and late_buffer == 3:

		ax.scatter(nucleus.centroid[2],nucleus.centroid[1],nucleus.centroid[0], s=10, c='green', marker='o')
		if text == True:
			ax.text(nucleus.centroid[2], nucleus.centroid[1], nucleus.centroid[0], f'{tp_late:03}_{nucleus.label:0{late_buffer}}')
		continue
	if f'{tp_late:03}_{nucleus.label:0{late_buffer}}' in only_these_nuclei:
		if f'{tp_late:03}_{nucleus.label:0{late_buffer}}' in highlight:
			x=np.append(x,nucleus.centroid[2])
			y=np.append(y,nucleus.centroid[1])
			z=np.append(z,nucleus.centroid[0])
			color=np.append(color,'lime')
			name=np.append(name, f'{tp_late:03}_{nucleus.label:0{late_buffer}}')
			
			txt = ax.text(nucleus.centroid[2], nucleus.centroid[1], nucleus.centroid[0], f'{tp_late:03}_{nucleus.label:0{late_buffer}}')
			txt.set_visible(False)
			early_txt.append(txt)
				#ax.text(nucleus.centroid[2], nucleus.centroid[1], nucleus.centroid[0], f'{tp_late:03}_{nucleus.label:0{late_buffer}}')
		else:
			x=np.append(x,nucleus.centroid[2])
			y=np.append(y,nucleus.centroid[1])
			z=np.append(z,nucleus.centroid[0])
			color=np.append(color,'blue')
			name=np.append(name,f'{tp_late:03}_{nucleus.label:0{late_buffer}}')
				
			txt = ax.text(nucleus.centroid[2], nucleus.centroid[1], nucleus.centroid[0], f'{tp_late:03}_{nucleus.label:0{late_buffer}}')
			txt.set_visible(False)
			early_txt.append(txt)
				#ax.text(nucleus.centroid[2], nucleus.centroid[1], nucleus.centroid[0], f'{tp_late:03}_{nucleus.label:0{late_buffer}}')
#print(x)
sc = ax.scatter(x,y,z, s=10, c=color, marker='o')


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
x_edges = np.array([])
y_edges = np.array([])
z_edges = np.array([])
edges_names = np.array([])
edge_lines = []
counter = 0
for edge in edges_to_show:
	ind_early = np.where(name == edge[0])	
	ind_late = np.where(name == edge[1])
	#print(np.array([x[ind_early][0], x[ind_late][0]]))
	if counter > 0:
		x_edges = np.vstack((x_edges, np.array([x[ind_early][0], x[ind_late][0]])))
		y_edges = np.vstack((y_edges, np.array([y[ind_early][0], y[ind_late][0]])))
		z_edges = np.vstack((z_edges, np.array([z[ind_early][0], z[ind_late][0]])))
		edges_names = np.vstack((edges_names, np.array([[name[ind_early].item(),name[ind_late].item()]])))
	else:
		x_edges = np.array([x[ind_early][0], x[ind_late][0]])
		y_edges = np.array([y[ind_early][0], y[ind_late][0]])
		z_edges = np.array([z[ind_early][0], z[ind_late][0]])
		edges_names = np.array([[name[ind_early].item(),name[ind_late].item()]])
		#print(edges_names.shape)
		
	
	edge_line, = ax.plot(x_edges[counter], y_edges[counter],z_edges[counter], color = 'white')
	edge_lines.append(edge_line)
	counter = counter + 1
#print(edges_names.shape)
#print(x_edges.shape)
#print(edges_names[4])



ax.set_facecolor('purple')
fig.patch.set_facecolor('black')

#print(len(edge_lines))
ax.set_xticks(range(0,(nuclear_seg_late.shape[2]), 100))
ax.set_yticks(range(0,(nuclear_seg_late.shape[1]), 30))
ax.set_zticks(range(0,(nuclear_seg_late.shape[0]), 100))



sc.set_picker(1)
last_text = None

visible = np.ones(len(x), dtype=bool)
edges_visible = np.ones(len(lines), dtype=bool)

idx_stack = []
clicks = []
click_idx = []


adding_edges = []
deleting_edges = []

def on_click(event):
	if event.button is MouseButton.RIGHT:
		mouse = np.array([event.x, event.y])
		dists = np.linalg.norm(projected - mouse, axis=1)
		idx = np.argmin(dists)
		#print('right click!')
		if dists[idx] < 5:
			visible[idx] = False
			sc._offsets3d = (x[visible],y[visible],z[visible])
			sc.set_color(color[visible])
			idx_stack.append(idx)
			fig.canvas.draw_idle()

	if event.button is MouseButton.LEFT:
		global edges_names
		mouse = np.array([event.x, event.y])
		dists = np.linalg.norm(projected - mouse, axis=1)
		idx1 = np.argmin(dists)
		#print('left click!')
		if dists[idx1] < 5:
			clicks.append(early_txt[idx1].get_text())
			click_idx.append(idx1)
			print('adding click')
		else:
			print('nothing selected')
			click_idx.clear()
			clicks.clear()
		if len(clicks) == 2:
			if int(clicks[0][:3]) > int(clicks[1][:3]):
				clicks.reverse()
				click_idx.reverse()
			print(clicks)
			if int(clicks[0][:3]) == int(clicks[1][:3]):
				clicks.clear()
				click_idx.clear()
				print('same')
				return()
			
			#print(int(clicks[0].get_text()[:3]))
			#print(int(clicks[1].get_text()[:3]))
			if np.any(np.all(edges_names == clicks, axis=1)): #Deleting Edge
				edge_idx = np.where(np.all(edges_names == clicks, axis=1))[0][0]
				edges_names[edge_idx] = np.nan
				edge_lines[edge_idx].remove()
				deleting_edges.append([early_txt[click_idx[0]].get_text(),early_txt[click_idx[1]].get_text()])
				print('deleted edge')
				
			else: #Adding Edge
				
				print(f'{early_txt[click_idx[0]].get_text()} position is {x[click_idx[0]]},{y[click_idx[0]]},{z[click_idx[0]]}')
				print(f'{early_txt[click_idx[1]].get_text()} position is {x[click_idx[1]]},{y[click_idx[1]]},{z[click_idx[1]]}')
				print('added edge')
				new_edge, = ax.plot([x[click_idx[0]],x[click_idx[1]]], [y[click_idx[0]],y[click_idx[1]]],[z[click_idx[0]],z[click_idx[1]]], color = 'white')
				edge_lines.append(new_edge)
				small = [early_txt[click_idx[0]].get_text(),early_txt[click_idx[1]].get_text()]
				print(small)
				edges_names = np.vstack((edges_names, small))
				adding_edges.append(small)

			clicks.clear()
			click_idx.clear()
			fig.canvas.draw_idle()




def on_key(event):
	if event.key == 'u':
		if len(idx_stack) > 0:
			my_idx = idx_stack.pop(len(idx_stack)-1)
			visible[my_idx] = True
			sc._offsets3d = (x[visible], y[visible],z[visible])	
			sc.set_color(color[visible])
			fig.canvas.draw_idle()

	if event.key == 'p':
		
		if len(adding_edges) == 0 or len(deleting_edges) == 0:
			print('nothing')
			print(f'adding: {adding_edges}')
			print(f'deleting: {deleting_edges}')
		elif len(adding_edges) == 0:
			print('only deleting')
			print(f'adding: {adding_edges}')
			print(f'deleting: {deleting_edges}')
		elif len(deleting_edges) == 0:
			print('only adding')
			print(f'adding: {adding_edges}')
			print(f'deleting: {deleting_edges}')
		else:


			
			while len(deleting_edges) != 0 and len(adding_edges) != 0:
				if len(np.vstack((adding_edges, deleting_edges))) == len(np.unique(np.vstack((adding_edges, deleting_edges)),axis=0)):
					break
				for edge in adding_edges:
					if edge in deleting_edges:
						adding_edges.remove(edge)
						deleting_edges.remove(edge)
			print(f'adding: {adding_edges}')
			print(f'deleting: {deleting_edges}')

		with open(f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/csvs/edges.csv', 'a', newline='') as f:
			writer = csv.writer(f)
			writer.writerows(adding_edges)
		print('saved adding edges')
		adding_edges.clear()
		new_csv = []
		with open(f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/csvs/edges.csv', 'r', newline='') as f:
			reader = csv.reader(f)
			for row in reader:
				if row not in deleting_edges:
					new_csv.append(row)
		with open(f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/csvs/edges.csv', 'w', newline='') as f:
			csv.writer(f).writerows(new_csv)
		deleting_edges.clear()
		print('saved deleting edges')
		new_graph = []
		idx_stack_new = [early_txt[x].get_text() for x in idx_stack]
		print(f'deleting: {idx_stack_new}')
		with open(f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/csvs/graph.csv', 'r', newline='') as f:
			reader = csv.reader(f)
			for row in reader:
				if row[0] not in idx_stack_new:
					new_graph.append(row)
		with open(f'/mnt/home/{user}/ceph/tracked_embryos/{movie}/csvs/graph.csv', 'w', newline='') as f:
			csv.writer(f).writerows(new_graph)
		print('saved deleting nuclei')

def hover(event):
	if event.inaxes!= ax:
		for i in early_txt:
			i.set_visible(False)
	mouse = np.array([event.x, event.y])
	dists = np.linalg.norm(projected - mouse, axis=1)
	idx = np.argmin(dists)
	
	if dists[idx] < 5:
		early_txt[idx].set_visible(True)
	else:
		for i in early_txt:
			i.set_visible(False)
	fig.canvas.draw_idle()

fig.canvas.mpl_connect("motion_notify_event", hover)
fig.canvas.mpl_connect("key_press_event", on_key)
fig.canvas.mpl_connect("button_press_event", on_click)
projected = None

#print(x_early)
def update_projection():
	global projected
	x2, y2, _ = proj3d.proj_transform(x,y,z,ax.get_proj())
	projected = ax.transData.transform(np.column_stack([x2,y2]))
	

fig.canvas.mpl_connect('draw_event', lambda event: update_projection())
plt.show()






	
	








