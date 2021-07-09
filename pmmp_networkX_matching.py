import pandas as pd
import matplotlib.pyplot as plt
import random
import numpy as np
import time

start = time.time()
#####Load data to a dataframe#####
df = pd.read_csv('data', delimiter = "\t")

#Get all different pairs of spectrum_key peptide_key values
samples = df[['spectrum_key','peptide_key']].drop_duplicates()
matchings_info = []
for index, sample in samples.iterrows():
    
    ######Get small example#####
    sample = df.loc[(df['spectrum_key'] == sample['spectrum_key']) & (df['peptide_key'] == sample['peptide_key'])]

    #####Create graph model######

    #number of modification sites
    n = len(sample['modification_site'].unique().tolist())
    #occurences of each modification type
    modifications_num_dict = {float(x.split(':')[0]) : int(x.split(':')[1]) for x in sample['n_mods'].iloc[1].split('_')}
    #number present of modifications
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
    V_names = ['site '+str(i) for i in A] + ['mod '+str(i)+' occ '+str(j) for i in range(len(keys)) for j in     range(modifications_num_dict[keys[i]])]
    #print('Vertices:')
    #for v in V_names:
        #print(v)

    #list of graph's edges
    E = []
    for index, row in sample.iterrows():
        net_edges = [(V_names[A.index(row['modification_site'])], V_names[len(A)+i], row['p_score']) for i in         np.where(np.array(D)==row['modification_mass'])[0].tolist()]
        E = E + net_edges
    #print('\nEdges:')
    #for e in E:
        #print(e)

    #####Use networkX library######

    from networkx.algorithms import bipartite
    import networkx as nx

    G = nx.Graph()
    # Add nodes with the node attribute "bipartite"
    G.add_nodes_from(V_names[:len(A)], bipartite=0)
    G.add_nodes_from(V_names[len(A):], bipartite=1)

    # Add edges only between nodes of opposite node sets
    G.add_weighted_edges_from(E)

    #Draw the grpah
    #X, Y = bipartite.sets(G)
    #pos = dict()
    #pos.update( (n, (1, i)) for i, n in enumerate(X) ) # put nodes from X at x=1
    #pos.update( (n, (2, i)) for i, n in enumerate(Y) ) # put nodes from Y at x=2

    #f, axs = plt.subplots(1,1,figsize=(24,20))
    #nx.draw(G, pos=pos, node_size=1000, with_labels=True, font_weight='bold')

    #Draw the weights of the edges
    #labels = nx.get_edge_attributes(G,'weight')
    #nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
    #plt.savefig('sample_plot.png')
    #plt.show()

    ####Find maximum weight matching#####
    matching = nx.algorithms.max_weight_matching(G)
    #print(matching)
    #Compute sum of weights of the matching
    matching_weights_sum = 0.0
    for m in matching:
        match_info = [sample['spectrum_key'].iloc[0], sample['peptide_key'].iloc[0]]
        for e in E:
            if m[0] in str(e) and m[1] in str(e):
                matching_weights_sum += e[2]
                if 'site' in m[0]:
                    match_info.append(int(m[0].split(' ')[1]))
                    match_info.append(D[V_names.index(m[1])-len(A)])
                else:
                    match_info.append(int(m[1].split(' ')[1]))
                    match_info.append(D[V_names.index(m[0])-len(A)])
        matchings_info.append(match_info)           

results_df = pd.DataFrame(matchings_info, columns = ['spectrum_key', 'peptide_key', 'modification_site', 'modification_mass'])
print(results_df)
print('\n')
results_df.to_csv('matchings.csv', index=False, sep ='\t')
end = time.time()
print('Total running time : ', end-start)

