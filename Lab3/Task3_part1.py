from linecache import getline
from platform import node
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def getDegPdf(G):
    pdf = nx.degree_histogram(G)
    pdf /= np.sum(pdf)
    return pdf

def getShortestPathLengths(G):
    shortest_paths = list(dict(nx.shortest_path_length(G)).values())
    shortest_paths = [list(sp.values()) for sp in shortest_paths]
    shortest_paths = np.hstack(shortest_paths)
    return shortest_paths

def getAverageShortestPathLength(G):
    sp = getShortestPathLengths(G)
    return np.average(sp)

def getAverageClusterCoeff(G):
    cl = list(nx.clustering(G).values())
    return np.average(cl)

def getLargestComponent(G):
    return G.subgraph((max(nx.connected_components(G), key=len))).copy()

def removeNodeFraction(node_indeces, frac):
    num_nodes = len(node_indeces)
    num_removed = int(frac*num_nodes)
    return node_indeces[num_removed:]
    

# Load graphs from files
graph_ecoli = nx.read_edgelist("ecoli_metabolic_net.txt", delimiter=",")
graph_kidney = nx.read_edgelist("human_kidney_protein_net.txt", delimiter=",")
graph_yeast = nx.read_edgelist("yeast_gene_net.txt", delimiter=",")

nodes_ecoli = list(nx.nodes(graph_ecoli))
nodes_kidney = list(nx.nodes(graph_kidney))
nodes_yeast = list(nx.nodes(graph_yeast))

diams_ecoli = []
diams_kidney = []
diams_yeast = []

# Add original graph diameter
diams_ecoli.append(max(getShortestPathLengths(graph_ecoli)))
diams_kidney.append(max(getShortestPathLengths(graph_kidney)))
diams_yeast.append(max(getShortestPathLengths(graph_yeast)))

fracs = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25]
samples = 20
for f in fracs[1:]:
    diam_ecoli = 0 
    diam_kidney = 0
    diam_yeast = 0
    for i in range(samples):
        # Sample fraction of nodes
        np.random.shuffle(nodes_ecoli)
        frac_nodes_ecoli = removeNodeFraction(nodes_ecoli, f)

        np.random.shuffle(nodes_kidney)
        frac_nodes_kidney = removeNodeFraction(nodes_kidney, f)
        
        np.random.shuffle(nodes_yeast)
        frac_nodes_yeast = removeNodeFraction(nodes_yeast, f)
        
        # Produce corresponding graphs
        frac_graph_ecoli = nx.subgraph(graph_ecoli, frac_nodes_ecoli)
        frac_graph_kidney = nx.subgraph(graph_kidney, frac_nodes_kidney)
        frac_graph_yeast = nx.subgraph(graph_yeast, frac_nodes_yeast)

        # Extract largest component of graph
        frac_graph_ecoli = getLargestComponent(frac_graph_ecoli)
        frac_graph_kidney = getLargestComponent(frac_graph_kidney)
        frac_graph_yeast = getLargestComponent(frac_graph_yeast)

        # Calculate diameters
        diam_ecoli += max(getShortestPathLengths(frac_graph_ecoli))
        diam_kidney += max(getShortestPathLengths(frac_graph_kidney))
        diam_yeast += max(getShortestPathLengths(frac_graph_yeast))

    diams_ecoli.append(diam_ecoli/samples)
    diams_kidney.append(diam_kidney/samples)
    diams_yeast.append(diam_yeast/samples)
        

print(diams_ecoli)
print(diams_kidney)
print(diams_yeast)

with open("data/rand_frac.txt", "w") as rand_frac_data:
    rand_frac_data.write(" ".join(str(f) for f in fracs) + "\n")
    rand_frac_data.write(" ".join(str(e) for e in diams_ecoli) + "\n")
    rand_frac_data.write(" ".join(str(k) for k in diams_kidney) + "\n")
    rand_frac_data.write(" ".join(str(y) for y in diams_yeast) + "\n")


    
