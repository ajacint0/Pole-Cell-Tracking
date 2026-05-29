def trim_graph(length, node, data, node_list, delete_list, cutoff):
	merges = []
	for edge in data:
		if edge[0] == node:
			merges.append(edge[1])

		if len(merges) == 2:
			break
	if node not in node_list:
		node_list.append(node)

	if len(merges) == 1:

		add_on1 = int(merges[0].split('_')[0]) - int(node.split('_')[0])

		node_list1, delete_list1 = trim_graph(length + add_on1, merges[0], data, node_list, delete_list, cutoff)

		node_list = node_list + node_list1
		delete_list = delete_list + delete_list1
			
	
	elif len(merges) == 2:

		add_on1 = int(merges[0].split('_')[0]) - int(node.split('_')[0])
		add_on2 = int(merges[1].split('_')[0]) - int(node.split('_')[0])

		node_list1, delete_list1 = trim_graph(length + add_on1, merges[0], data, node_list, delete_list, cutoff)
		node_list2, delete_list2 = trim_graph(length + add_on2, merges[1], data, node_list, delete_list, cutoff)

		node_list = node_list + node_list1 + node_list2
		delete_list = delete_list + delete_list1 + delete_list2
		
	
	else:
		
		if length < cutoff:
			delete_list.append(True)
			
		else:
			delete_list.append(False)

	return node_list, delete_list
