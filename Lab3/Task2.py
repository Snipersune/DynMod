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


# Load graphs from files
graph_ecoli = nx.read_edgelist("ecoli_metabolic_net.txt", delimiter=",")
graph_kidney = nx.read_edgelist("human_kidney_protein_net.txt", delimiter=",")
graph_yeast = nx.read_edgelist("yeast_gene_net.txt", delimiter=",")

# Get degree probability density functions
deg_pdf_ecoli = getDegPdf(graph_ecoli)
deg_pdf_kidney = getDegPdf(graph_kidney)
deg_pdf_yeast = getDegPdf(graph_yeast)


# Degrees in respective distribution
degs_ecoli = np.arange(0, len(deg_pdf_ecoli))
degs_kidney = np.arange(0, len(deg_pdf_kidney))
degs_yeast = np.arange(0, len(deg_pdf_yeast))

avg_deg_ecoli = np.dot(deg_pdf_ecoli, degs_ecoli)
avg_deg_kidney = np.dot(deg_pdf_kidney, degs_kidney)
avg_deg_yeast = np.dot(deg_pdf_yeast, degs_yeast)

num_nodes_ecoli = nx.number_of_nodes(graph_ecoli)
num_nodes_kidney = nx.number_of_nodes(graph_kidney)
num_nodes_yeast = nx.number_of_nodes(graph_yeast)

print(avg_deg_ecoli)
print(avg_deg_kidney)
print(avg_deg_yeast)

p_ecoli = avg_deg_ecoli/(num_nodes_ecoli-1)
p_kidney = avg_deg_kidney/(num_nodes_kidney-1)
p_yeast = avg_deg_yeast/(num_nodes_yeast-1)

# To store average cluster coeffs and average shortest 
# path lengths of erdos networks
avg_cl_erdos_ecoli = []
avg_sp_erdos_ecoli = []

avg_cl_erdos_kidney = []
avg_sp_erdos_kidney = []

avg_cl_erdos_yeast = []
avg_sp_erdos_yeast = []

# Generate Erdös-Reyni networks
for i in range(80):
    er_ecoli = nx.erdos_renyi_graph(num_nodes_ecoli, p_ecoli)
    er_kidney = nx.erdos_renyi_graph(num_nodes_kidney, p_kidney)
    er_yeast = nx.erdos_renyi_graph(num_nodes_yeast, p_yeast)
    
    avg_cl_erdos_ecoli.append(getAverageClusterCoeff(er_ecoli))
    avg_sp_erdos_ecoli.append(getAverageShortestPathLength(er_ecoli))
    
    avg_cl_erdos_kidney.append(getAverageClusterCoeff(er_kidney))
    avg_sp_erdos_kidney.append(getAverageShortestPathLength(er_kidney))

    avg_cl_erdos_yeast.append(getAverageClusterCoeff(er_yeast))
    avg_sp_erdos_yeast.append(getAverageShortestPathLength(er_yeast))
    


# Averages of actual networks
avg_cl_ecoli = getAverageClusterCoeff(graph_ecoli)
avg_sp_ecoli = getAverageShortestPathLength(graph_ecoli)

avg_cl_kidney = getAverageClusterCoeff(graph_kidney)
avg_sp_kidney = getAverageShortestPathLength(graph_kidney)

avg_cl_yeast = getAverageClusterCoeff(graph_yeast)
avg_sp_yeast = getAverageShortestPathLength(graph_yeast)



plt.figure()
plt.vlines(avg_cl_ecoli, ymin=0, ymax=10, colors="red")
plt.hist(avg_cl_erdos_ecoli, bins=10, label="ecoli")
plt.legend()

plt.figure()
plt.vlines(avg_cl_kidney, ymin=0, ymax=10, colors="red")
plt.hist(avg_cl_erdos_kidney, bins=10, label="kidney")
plt.legend()

plt.figure()
plt.vlines(avg_cl_yeast, ymin=0, ymax=10, colors="red")
plt.hist(avg_cl_erdos_yeast, bins=10, label="yeast")
plt.legend()

plt.show()

plt.figure()
plt.vlines(avg_sp_ecoli, ymin=0, ymax=10, colors="red")
plt.hist(avg_sp_erdos_ecoli, bins=10, label="ecoli")
plt.legend()

plt.figure()
plt.vlines(avg_sp_kidney, ymin=0, ymax=10, colors="red")
plt.hist(avg_sp_erdos_kidney, bins=10, label="kidney")
plt.legend()

plt.figure()
plt.vlines(avg_sp_yeast, ymin=0, ymax=10, colors="red")
plt.hist(avg_sp_erdos_yeast, bins=10, label="yeast")
plt.legend()

plt.show()
