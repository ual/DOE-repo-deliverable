# This function converts uses the bayarea_tertiary graphml file to create network files that igraph
# can interpret. This involvels adding labes from 1 up N to all the nodes, that will be used in place
# of the osmids. The conversion from igraph id back to osmid can be found in the igraph_node_labels.csv file


import networkx as nx
import csv
import numpy as np

def main():
    #Loading the graph from file into variable g
    #g = igraph.Graph.Read_GraphML("bayarea_tertiary_simplified.graphml")
    nxgraph = nx.read_graphml("bayarea_tertiary_simplified.graphml")

    # Add labels from 1 up N to the nodes
    label_dict = dict(zip(nxgraph.nodes.keys(),range(1,len(nxgraph.node.keys())+1)))
    nx.set_node_attributes(nxgraph, label_dict, name='igraph_label')

    #First write to csv the new labels to be used later

    with open("igraph_node_labels.csv", "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['osmid','igraph_id'])
        for key, value in label_dict.items():
            row = [int(key), value]
            writer.writerow(row)


    #Loading the demand
    demand_Location = "pm_peak_nodes.csv"
    demand = np.loadtxt(demand_Location, delimiter=',', skiprows=1)

    # Change demand to use igraph labels as well
    csv_file = open("igraph_demand.csv","w", newline="")
    writer = csv.writer(csv_file)
    writer.writerow(['orig','dest','trips'])
    for i in range(demand.shape[0]):
        row = list(demand[i,:])
        row[0] = int(nxgraph.node[str(int(demand[i,0]))]['igraph_label'])
        row[1] = int(nxgraph.node[str(int(demand[i,1]))]['igraph_label'])
        writer.writerow(row)

    #Load BPR functions
    bpr_Location = "bay_area_bpr.csv"
    bpr = np.loadtxt(bpr_Location, delimiter=',', skiprows=1)


    #Rewrite the bpr file
    csv_file = open("igraph_bpr_coefficients.csv","w", newline="")
    writer = csv.writer(csv_file)
    writer.writerow(['LINK','A','B','a0','a1','a2','a3','a4'])
    for i in range(bpr.shape[0]):
        row = list(bpr[i,:])
        row[0] = int(row[0])
        row[1] = int(nxgraph.node[str(int(bpr[i,1]))]['igraph_label'])
        row[2] = int(nxgraph.node[str(int(bpr[i,2]))]['igraph_label'])
        print(row)
        writer.writerow(row)

if __name__ == "__main__": main()