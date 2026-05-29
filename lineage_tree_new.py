import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import tifffile as tif
import csv
from trim_graph import trim_graph



limit_graph = True
asymmetric = ['048_000014','048_000623', '048_000723', '049_000544']
cutoff = 20
G = nx.Graph()
cust_pos = {}

arr = np.load('/mnt/home/ajacinto/ceph/tracked_embryos/2025-12-17_160948/important_nuclei/tp_35_for_sure.npy')
arr_maybe = np.load('/mnt/home/ajacinto/ceph/tracked_embryos/2025-12-17_160948/important_nuclei/tp_35_maybe.npy')

#arr = []
#MAYBE RETURN 0 IF THE LINEAGE IS TO BE DELETED, THEN DELETE THE NODES THAT HAVE BEEN SAVED IN A LIST

def rec(depth, node, csvfile, line_count):
	if node == '083_001115':
		print('hi')
	merges = []

	if node in cust_pos:
		pass
	else:	
		
		cust_pos.update({node: (int(node.split('_')[0]) / 10,line_count)})

	for edge in data:
		if edge[0] == node:
			merges.append(edge[1])

		if len(merges) == 2:
			break
	if len(merges) == 1:
		rec(depth, merges[0], csvfile, line_count)

	elif len(merges) == 2:
		depth = depth + 1
		up_line_count = 200/ (2**depth)
		down_line_count = 200 / (2**depth) * -1
		rec(depth, merges[0], csvfile, line_count + up_line_count)
		rec(depth, merges[1], csvfile, line_count + down_line_count)
	else:
		return
		
		




with open('/mnt/home/ajacinto/ceph/tracked_embryos/2025-12-17_160948/csvs/volume_graph.csv', newline='') as file:
	csvfile = csv.reader(file)
	for lines in csvfile:
		G.add_node(lines[0], intensity=lines[1], nucleus=lines[0].split('_')[0])
		
with open('/mnt/home/ajacinto/ceph/tracked_embryos/2025-12-17_160948/csvs/volume_edges.csv', newline='') as file:
	csvfile = csv.reader(file)
	data = list(csvfile)

	for lines in data:
		G.add_edge(lines[0],lines[1])

line_count = 0


for i in range(15, 114):
	with open('/mnt/home/ajacinto/ceph/tracked_embryos/2025-12-17_160948/csvs/volume_edges.csv', newline='') as file:
		csvfile = csv.reader(file)
		
		for lines in csvfile:
			if int(lines[0].split('_')[0]) != i:
				continue
			
			if lines[0] not in cust_pos:
				
				
				rec(0, lines[0], data, line_count)
				line_count = line_count + 800


if limit_graph == True:
	for i in range(15, 114):
		
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
				if False not in delete_list:
					for vertex in node_list:
						try:
							G.remove_node(vertex)
						except:
							continue
					
				#print(node, delete_list)
		



fig = plt.figure()
nx.draw_networkx_nodes(G, pos=cust_pos, node_size = 20)
nx.draw_networkx_edges(G, pos=cust_pos, width=2)
nx.draw_networkx_labels(G, pos=cust_pos, font_size=9, font_family='sans-serif', verticalalignment='top')
#plt.set_facecolor('purple')

for node in arr:
	nx.draw_networkx_nodes(G, pos={node: cust_pos[node]}, nodelist = [node], node_size = 20, node_color = '#FFFF00')

for node in arr_maybe:
	nx.draw_networkx_nodes(G, pos={node: cust_pos[node]}, nodelist = [node], node_size = 20, node_color = '#FFA500')



for node in asymmetric:
	nx.draw_networkx_nodes(G, pos={node: cust_pos[node]}, nodelist = [node], node_size = 20, node_color = '#FF0000')



fig.patch.set_facecolor('purple')
plt.axis('off')
plt.show()
