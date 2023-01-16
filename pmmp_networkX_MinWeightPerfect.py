import pandas as pd
import matplotlib.pyplot as plt
import random
import numpy as np
import time

class LessSitesException(Exception):
        "Raised when the possible amino acid chain sites are less than the modifications."
        pass
        
class DisconnectedGraphException(Exception):
        "Raised when there is at least one modification without an edge to a possible site in the amino acid chain."
        pass

start = time.time()
#####Load data to a dataframe#####
df = pd.read_csv('./datasets/Phospho/data', delimiter = "\t")

#Get all different pairs of spectrum_key peptide_key values
samples = df[['spectrum_key','peptide_key']].drop_duplicates()
#df.loc[df['p_score'] == 0.0, 'p_score'] = -1
matchings_info = []
for index, sample in samples.iterrows():
    
    ######Get small example#####
    sample = df.loc[(df['spectrum_key'] == sample['spectrum_key']) & (df['peptide_key'] == sample['peptide_key'])]
    #####Create graph model######

    #number of modification sites
    n = len(sample['modification_site'].unique().tolist())
    #occurences of each modification type
    modifications_num_dict = {np.float16(x.split(':')[0]) : int(x.split(':')[1]) for x in sample['n_mods'].iloc[1].split('_')}
    #total number of modifications
    k = sum(modifications_num_dict.values())
    #list of vertices representing modification sites
    A = sample['modification_site'].unique().tolist()
    #list of vertices representing modifications
    D = [key for key in modifications_num_dict for i in range(modifications_num_dict[key])]
 
    #list of graph's vertices
    #vertex names: 
    #      'site 10' denotes the 10-th site of the amino acid chain
    #      'mod 1 occ 2' denote the second occurence of the modification of the type 1 (e.g. Oxidation)
    V = A + D
    keys = [key for key in modifications_num_dict]
    V_names = ['site '+str(i) for i in A] + ['mod '+str(i)+' occ '+str(j+1) for i in range(len(keys)) for j in     range(modifications_num_dict[keys[i]])]
    
    #Add dummy vertices to make a perfect matching possible
    if len(D) < len(A):
    	dummy_vertices = ['dummy_mod_' + str(i) for i in range(len(A)-len(D))]
    	V_names = V_names + dummy_vertices
    	
    #print('Vertices: ', len(V_names))
    #for v in V_names:
    #    print(v)

    #list of graph's edges
    E = []
    for index, row in sample.iterrows():
        net_edges = [[V_names[A.index(row['modification_site'])], V_names[len(A)+i], row['p_score']] for i in         np.where(np.array(D)==np.float16(row['modification_mass']))[0].tolist()]
        E = E + net_edges
        
    #print('Number of edges: ', len(E))
    #for e in E:
    #    print(e)
        
    #Normalize edge weights
    weights = np.array([e[2] for e in E])
    weights[weights == 0.0] = 0.01
    normalized_weights = weights / np.linalg.norm(weights)
    max_weight = max(normalized_weights) + 0.5
    reverted_weights = np.array([max(normalized_weights)-w for w in normalized_weights])
    reverted_weights[reverted_weights == 0.0] = 0.01

    for i in range(len(E)):
    	edge = E[i]
    	edge[2] = reverted_weights[i]
 
    #Make graph complete
    # add edges from/to dummy vertices
    if len(D) < len(A):
    	for a_vertex_index in range(len(A)):
    		for dummy_vertex_index in range(len(A) + len(D), len(V_names)):
    			dummy_edge = [V_names[a_vertex_index], V_names[dummy_vertex_index], max_weight]
    			E.append(dummy_edge)
    
    #Add missing edges from other vertices
    for v1_index in range(len(A)):
    	vertices_connected_to_v1 = [e[1] for e in E if (e[0] == V_names[v1_index])]
    	for v2_index in range(len(A)+1, len(A) + len(D)):
    		if V_names[v2_index] in vertices_connected_to_v1:
    			continue
    		else:
    			edge = [V_names[v1_index], V_names[v2_index], max_weight]

    #check if the graph does not satisfy requirements
    try:
        if len(A) < len(D):
            raise LessSitesException
	
    except LessSitesException:
        print("Exception occurred: less sites than modifications")
        continue

    #####Use networkX library######

    from networkx.algorithms import bipartite
    import networkx as nx

    G = nx.Graph()
    # Add nodes with the node attribute "bipartite"
    G.add_nodes_from(V_names[:len(A)], bipartite=0)
    G.add_nodes_from(V_names[len(A):], bipartite=1)

    # Add edges only between nodes of opposite node sets
    G.add_weighted_edges_from(E)
    
    #check if the graph does not satisfy requirements
    try:
        connected = nx.is_connected(G)
        if connected == False:
            raise DisconnectedGraphException
	
    except DisconnectedGraphException:
        print("Exception occurred: there are modifications without a possible amino acid chain site.")
        continue

    

    ####Find maximum weight matching#####
    matching = nx.algorithms.bipartite.matching.minimum_weight_full_matching(G)
    #Compute sum of weights of the matching
    matching_weights_sum = 0.0
    matched = []
    for key in matching:
        vertex1 = key
        vertex2 = matching[key]
        if vertex1 in matched or vertex2 in matched:
            continue
        match_info = [sample['spectrum_key'].iloc[0], sample['peptide_key'].iloc[0]]
        flag_dummy_vertex = False
        for e in E:
            if vertex1 in e and vertex2 in e:
                if 'dummy' in vertex1 or 'dummy' in vertex2:
                    flag_dummy_vertex = True
                    break
                matching_weights_sum += e[2]
                if 'site' in vertex1:
                    matched.append(vertex1)
                    matched.append(vertex2)
                    match_info.append(int(vertex1.split(' ')[1]))
                    match_info.append(D[V_names.index(vertex2)-len(A)])
                else:
                    match_info.append(int(vertex2.split(' ')[1]))
                    match_info.append(D[V_names.index(vertex1)-len(A)])
        if flag_dummy_vertex == True:
        	continue
        print(match_info)
        matchings_info.append(match_info)
        
    #Draw the grpah
    #X, Y = bipartite.sets(G)
    #pos = dict()
    #pos.update( (n, (1, i)) for i, n in enumerate(X) ) # put nodes from X at x=1
    #pos.update( (n, (2, i)) for i, n in enumerate(Y) ) # put nodes from Y at x=2

    #f, axs = plt.subplots(1,1,figsize=(24,20))
    #nx.draw(G, pos=pos, node_size=1000, with_labels=True, font_weight='bold')

    #plt.show()

results_df = pd.DataFrame(matchings_info, columns = ['spectrum_key', 'peptide_key', 'modification_site', 'modification_mass'])
print(results_df)
print('\n')
results_df.to_csv('./matchings_phospho.csv', index=False, sep ='\t')
end = time.time()
print('Total running time : ', end-start)

