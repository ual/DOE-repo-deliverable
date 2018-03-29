import networkx as nx
import csv
import numpy as np

def main():
    #Loading the graph from file into variable g
    #g = igraph.Graph.Read_GraphML("bayarea_tertiary_simplified.graphml")
    nxgraph = nx.read_graphml("bayarea_tertiary_simplified.graphml")

    demand = np.loadtxt('no_paths.csv', delimiter=',', skiprows=1)
    net_to_igraph = np.loadtxt('bay_area_ter_igraph_node_labels.csv', delimiter=',', skiprows=1)
    conversion_dic = dict(zip(net_to_igraph[:,1],net_to_igraph[:,0]))

    for i in range(0, demand.shape[0]):
        origin = demand[i,0]
        o_id = int(conversion_dic[origin])
        destination = demand[i,1]
        d_id = int(conversion_dic[origin])
        path = nx.shortest_path(nxgraph, source=str(o_id), target=str(d_id), weight=None)
        print (path)

if __name__ == "__main__": main()