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

# Get graphs from files
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

# Plot the results
plt.figure()
plt.scatter(degs_ecoli, deg_pdf_ecoli, label="ecoli")
plt.scatter(degs_kidney, deg_pdf_kidney, label="kidney")
plt.scatter(degs_yeast, deg_pdf_yeast, label="yeast")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Degree $k$")
plt.ylabel("$P(k)$")
plt.legend()
plt.show()


cl_ecoli = nx.clustering(graph_ecoli).values()
cl_kidney = nx.clustering(graph_kidney).values()
cl_yeast = nx.clustering(graph_yeast).values()

plt.figure()
plt.hist(cl_ecoli, bins=10, label="ecoli")
plt.legend()

plt.figure()
plt.hist(cl_kidney, bins=10, label="kidney")
plt.legend()

plt.figure()
plt.hist(cl_yeast, bins=10, label="yeast")
plt.legend()

plt.show()

sp_ecoli = getShortestPathLengths(graph_ecoli)
sp_kidney = getShortestPathLengths(graph_kidney)
sp_yeast = getShortestPathLengths(graph_yeast)

plt.figure()
plt.hist(sp_ecoli, bins=np.arange(0,10)-0.5, label="ecoli")
plt.legend()

plt.figure()
plt.hist(sp_kidney, bins=np.arange(0,10)-0.5, label="kidney")
plt.legend()

plt.figure()
plt.hist(sp_yeast, bins=np.arange(0,10)-0.5, label="yeast")
plt.legend()

plt.show()
