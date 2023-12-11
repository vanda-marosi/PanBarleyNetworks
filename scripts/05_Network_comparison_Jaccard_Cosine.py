#!/usr/bin/env python
# coding: utf-8

## Calculate Jaccard distance & Cosine similarity across co-expression networks
import pandas as pd
import numpy as np
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import csv
import os
import re
import glob
import pickle
import collections
import math

# Import WGCNA module tables
os.getcwd()
fnames = glob.glob('*_ortmodnames.csv')
# Initialize an empty dictionary for the merged results
merged_dict = {}

for fname in fnames:
    df = pd.read_csv(fname, delimiter = ",")
    df = df[['module_name', 'ID']]
    print(df.shape)
    mod_ort = df.groupby('module_name')['ID'].apply(set).to_dict()
    print(len(mod_ort))
    merged_dict = {**merged_dict, **mod_ort} # Merge the current dictionary into the merged dictionary
    
# Print the merged dictionary
print(len(merged_dict)) 

###########################
### Save nodelist for next step networkx analysis
nodes = list(merged_dict.keys())
print(len(nodes))
print(nodes[0:3])
# Write to pickle file
with open('all_nodes.pkl', 'wb') as f:
    pickle.dump(nodes, f)

### Save modules' size-table for next step networkx analysis
# Import WGCNA module tables
os.getcwd()
fnames = glob.glob('*_ortmodnames.csv')
node_size = {}
for fname in fnames:
    df = pd.read_csv(fname, delimiter = ",")
    df = df[['module_name', 'size']]
    mod_ort = df.groupby('module_name')['size'].first().to_dict()
    node_size = {**node_size, **mod_ort}
# Write to pickle file
with open('nodes_size.pkl', 'wb') as f:
    pickle.dump(node_size, f)
############################

#####################
### Calculate Jaccard distances
#####################
# Function for Jaccard similarity
def jaccard_similarity(set1, set2):
    """Compute the Jaccard similarity between two sets."""
    u = set1.union(set2)
    i = set1.intersection(set2)
    return len(i) / len(u)

# Compute the Jaccard coefficients for all unique pairs of modules.
calc_jaccard = {}
jaccard_list = []
for module1, module2 in itertools.product(merged_dict, repeat = 2):
        key = f"{module1}-{module2}"
        jaccard = jaccard_similarity(merged_dict[module1], merged_dict[module2])
        calc_jaccard[key] = round(jaccard, 4)
        jaccard_list.append(round(jaccard, 4))

# Print an intersected example    
print(dict(list(calc_jaccard.items())[0:5]))
# Number of intersected modules (len(d1)*len(d2))
print(len(calc_jaccard))
print(jaccard_list[0:10])
# Save Jaccard list
with open("jaccard_dict.bin", "wb") as output:
    pickle.dump(calc_jaccard, output)

### Plot Jaccard distance values on heatmap
############################
# Generate a symmetric matrix and dataframe
chunked_list = list()
chunk_size = len(merged_dict) # module length
for i in range(0, len(jaccard_list), chunk_size):
    chunked_list.append(jaccard_list[i:i+chunk_size])
print(chunked_list[0:2])
print(len(chunked_list))
# Get module name pair combinations into string list to use as colnames
module_names = list()
for key in merged_dict.keys():
    module_names.append(key)
print(len(module_names))
# Combine above object into df
jaccard_df = pd.DataFrame(chunked_list, columns =module_names)
jaccard_df.index = module_names
print(jaccard_df.shape)
print(jaccard_df[0:1])
# Clustered heatmap of Jaccard values
sns.set_context('talk', font_scale=2)
jaccard_clust = sns.clustermap(jaccard_df, cmap="mako_r", xticklabels=False, yticklabels=False, 
                               method="average", figsize=(11, 10), dendrogram_ratio=(.05, .05))
# Adjust colorbar position and size
jaccard_clust.cax.set_position([0.95, .2, .015, .3])  # x, y, width, height
jaccard_clust.cax.set_title('Jaccard-index', fontsize=12)  # Adjust colorbar title size
jaccard_clust.cax.tick_params(labelsize=8)  # Adjust colorbar tick label size
plt.tight_layout()  # Adjusts layout so labels don't overlap with other elements
plt.savefig('Jaccard_heatmap.pdf', format='pdf')
plt.show()
# Histogram of Jaccard values
bins = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
plt.hist(jaccard_list, bins=bins, edgecolor="k", alpha=0.7)
plt.yticks(fontsize=12)  
plt.xticks(bins, fontsize=12, rotation=45)  # Adjust fontsize and rotation for clarity
plt.xlabel("Jaccard Coefficient", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
plt.tight_layout()  # Adjusts layout so labels don't overlap with other elements
plt.savefig('Jaccard_histogram.pdf', format='pdf')
plt.show()

####################
## Calculate Cosine-similarity
####################
def word2vec(word):
    from collections import Counter
    from math import sqrt
    # count the characters in word
    cw = Counter(word)
    # precomputes a set of the different characters
    sw = set(cw)
    # precomputes the "length" of the word vector
    lw = sqrt(sum(c*c for c in cw.values()))
    # return a tuple
    return cw, sw, lw

def cosdis(v1, v2):
    # which characters are common to the two words?
    common = v1[1].intersection(v2[1])
    # by definition of cosine distance we have
    return sum(v1[0][ch]*v2[0][ch] for ch in common)/v1[2]/v2[2]

def precompute_vectors(merged_dict):
    """Precompute word vectors for all modules."""
    vectors = {}
    for module, ids in merged_dict.items():
        vectors[module] = word2vec(list(ids))
    return vectors

# Precompute the vectors for each module
precomputed_vectors = precompute_vectors(merged_dict)

# Adjusted Cosine distance calculation using precomputed vectors
calc_cosine = {}
cosine_list = []
for module1, module2 in itertools.product(merged_dict, repeat=2):
    key = f"{module1}-{module2}"
    cosine = cosdis(precomputed_vectors[module1], precomputed_vectors[module2])
    calc_cosine[key] = round(cosine, 4)
    cosine_list.append(round(cosine, 4))
# Print an example    
dict(list(calc_cosine.items())[0:20]), len(calc_cosine), cosine_list[0:10]

# Save Cosine list
with open("cosine_dict.bin", "wb") as output:
    pickle.dump(calc_cosine, output)

### Plot Cosine similarity values on heatmap
############################
# Split list into chunks using For Loop
chunked_list = list()
chunk_size = len(merged_dict) # module length
for i in range(0, len(cosine_list), chunk_size):
    chunked_list.append(cosine_list[i:i+chunk_size])
print(chunked_list[0:2])
print(len(chunked_list))

# Get module name pair combinations into str list to use as column names
module_names = list()
for key in merged_dict.keys():
    module_names.append(key)
print(module_names)

# Combine above object into df
cosine_df = pd.DataFrame(chunked_list, columns =module_names)
cosine_df.index = module_names
print(cosine_df.shape)
print(cosine_df[0:1])
# Save table
cosine_df.to_csv("cosine_df.csv", index=True)

# Clustered heatmap of Cosine values
sns.set_context('talk', font_scale=2)
jaccard_clust = sns.clustermap(cosine_df, cmap="mako_r", xticklabels=False, yticklabels=False, 
                               method="average", figsize=(11, 10), dendrogram_ratio=(.05, .05))

# Add lines to separate every 20 items
ax = jaccard_clust.ax_heatmap
for i in range(0, 800, 20):
    ax.axhline(i, color='white', lw=1)
    ax.axvline(i, color='white', lw=1)

# Adjust colorbar position and size
jaccard_clust.cax.set_position([0.95, .2, .015, .3])  # x, y, width, height
jaccard_clust.cax.set_title('Cosine', fontsize=12)  # Adjust colorbar title size
jaccard_clust.cax.tick_params(labelsize=8)  # Adjust colorbar tick label size
plt.tight_layout()  # Adjusts layout so labels don't overlap with other elements
plt.savefig('Cosine_heatmap.pdf', format='pdf')
plt.show()

# Histogram
bins = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
plt.hist(cosine_list, bins=bins, edgecolor="k", alpha=0.7)
plt.yticks(fontsize=12)  
plt.xticks(bins, fontsize=12, rotation=45)  # Adjust fontsize and rotation for clarity
plt.xlabel("Cosine similarity", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
plt.tight_layout()  # Adjusts layout so labels don't overlap with other elements
plt.savefig('Cosine_histogram.pdf', format='pdf')
plt.show()
