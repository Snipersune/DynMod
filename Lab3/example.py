import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# Read the network file as an edgelist, as specified in instructions
# Note that the delimiter is ","
# https://networkx.org/documentation/stable/reference/readwrite/generated/networkx.readwrite.edgelist.read_edgelist.html
G = nx.read_edgelist("human_kidney_protein_net.txt", delimiter=",")

# You can get the degree from every node as a dict.
# We use list comprehension to pick out just the degrees, ignoring the keys, which are just the node names.
# https://networkx.org/documentation/stable/reference/classes/generated/networkx.Graph.degree.html
degrees = [d for _, d in G.degree]

# We can compute the degree distribution by picking out the unique degrees and counting
udegrees = np.unique(degrees)
Pd = np.array([sum(degrees == ud) for ud in udegrees])
Pd = Pd/sum(Pd)

# Plot the results
plt.scatter(udegrees, Pd)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Degree $k$")
plt.ylabel("$P(k)$")
plt.show()