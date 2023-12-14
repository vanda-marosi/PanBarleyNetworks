#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
np.random.seed(108)
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
from netgraph import Graph
from collections import defaultdict
import pickle

# set matplotlib globals
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype'] = 42

os.getcwd()
with open("cosine_dict.bin", "rb") as data:
    cosine = pickle.load(data)
with open("all_nodes.pkl", "rb") as data:
    nodes = pickle.load(data)
with open("nodes_size.pkl", "rb") as data:
    nodes_size = pickle.load(data)

print("Package and data imports were successfull.")

# Separate edges based on weight values
edges = []
for key, item in cosine.items():
    if 0.07 <= item < 1.0: 
        nodes_key = key.split("-")
        wnodes = (nodes_key[0], nodes_key[1], item)
        edges.append(wnodes)
    else:
        pass
# Build graph
G = nx.Graph()
G.add_nodes_from(nodes)
G.add_weighted_edges_from(edges)
nx.set_node_attributes(G, nodes_size, 'size')
edge_widths = [G[u][v]['weight'] for u, v in G.edges()]

# Remove isolated nodes (nodes without any edges)
isolated_nodes = list(nx.isolates(G))
G.remove_nodes_from(isolated_nodes)
print(nx.number_connected_components(G))

# Compute communities
communities = list(nx.algorithms.community.louvain_communities(G, seed = 108, resolution = 1))
communities = sorted(communities, key=len, reverse=True)
# Assign each node to a community in a dictionary
partition = {}
for i, community in enumerate(communities):
    #print(i, community)
    for node in community:
        partition[node] = i
print(len(communities))
##############
## Save results in table
##############
meta_mod_table = []
for meta, modlist in partition.items():
    for mod in modlist:
        x = (mod, meta)
        meta_mod_table.append(x)
mod_meta_df = pd.DataFrame(meta_mod_table, columns =['Module', 'Community'])
print(mod_meta_df[0:5]) 
os.getcwd()
mod_meta_df.to_csv('Module_to_metamodule_louvain.csv', sep=',', index=False, header=False)
##############
print("Network is developed and communities are detected. Now figure development starts.")
##############
## Draw figure
##############
# Function to rescale the size values in the dictionary to a new range to fit for plotting
def rescale_dict_values(dictionary, new_min, new_max):
    old_min = min(dictionary.values())
    old_max = max(dictionary.values())
    # Rescaling formula: new_value = ((old_value - old_min) / (old_max - old_min)) * (new_max - new_min) + new_min
    scaled_dict = {
        key: ((value - old_min) / (old_max - old_min)) * (new_max - new_min) + new_min
        for key, value in dictionary.items()
    }
    return scaled_dict

# Define new range
new_min = 0.5
new_max = 3.0

# Rescale the dictionary values
rescaled_size = rescale_dict_values(nodes_size, new_min, new_max)

# color assignment
community_to_color = {
    0 : "#c7eae5",
    1 : "#5ab4ac",
    2 : "#d8b365",
    3 : "#01665e",
    4 : "#f6e8c3",
    5 : "#8c510a"
}
node_color = {node: community_to_color[community_id] \
              for node, community_id in partition.items()}

# Convert Community list of lists into dictionary
def create_community_dict(communities_list):
    # Define some properties to cycle through for the communities
    colors = ['#c7eae5', '#5ab4ac', '#d8b365', '#01665e', '#f6e8c3', '#8c510a']
    shapes = ['o', 'o', 'o', 'o', 'o', 'o']
    sizes = [10, 10, 10, 10, 10, 10]
    communities = {}
    for i, community_nodes in enumerate(communities_list):
        community_name = f'C{i+1}'
        communities[community_name] = {
            'nodes': community_nodes,
            'color': colors[i % len(colors)],
            'shape': shapes[i % len(shapes)],
            'size': sizes[i % len(sizes)]
        }
    return communities
communities_dict = create_community_dict(communities)
print(dict(list(communities_dict.items())[4:5]))

# Create proxy artists for legend handles
community_proxy_artists = []
for community, properties in communities_dict.items():
    proxy = plt.Line2D(
        [], [],
        linestyle='None',
        color=properties['color'],
        marker=properties['shape'],
        markersize=properties['size'],
        label=community
    )
    community_proxy_artists.append(proxy)
# Plot
fig, ax = plt.subplots(figsize=(10, 8))
Graph(G,
      node_color=node_color, # indicates the community each belongs to  
      node_edge_width=0.07,     # no black border around nodes 
      node_size = rescaled_size,       # radius, if dict, maps each node to an individual size, scalable
      node_alpha = 0.7,
      edge_width=0.2,        # use thin edges, as they carry no information in this visualisation
      edge_alpha=0.5,        # low edge alpha values accentuates bundles as they appear darker than single edges
      node_layout='community', node_layout_kwargs=dict(node_to_community=partition),
      edge_layout='bundled', edge_layout_kwargs=dict(k=2000), # this is where bundling is made possible
      ax=ax,
)

# Create legend for communities
community_legend = ax.legend(handles=community_proxy_artists, loc='upper right', title='Communities')
ax.add_artist(community_legend)
# Save
plt.savefig("Communities.png", dpi=300)
fig.savefig('Communities.svg', format='svg')
plt.savefig("Communities.pdf", format="pdf")

print("Network drawing was generated and saved successfully.")
