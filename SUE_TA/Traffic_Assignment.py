import igraph
import csv
import numpy as np
import timeit

# Function to load the graph from file
def load_graph(path_to_graph_file):
    g = igraph.Graph.Read_GraphML(path_to_graph_file)
    return g

def construct_igraph(graph):
    # 'vertices' contains the range of the vertices' indices in the graph
    x = int(np.min(graph[:,1:3]))
    y = int(np.max(graph[:,1:3]))+1
    vertices = range(x, y)
    # 'edges' is a list of the edges (to_id, from_id) in the graph
    edges = graph[:,1:3].astype(int).tolist()
    g = igraph.Graph(vertex_attrs={"label":vertices}, edges=edges, directed=True)
    g.es["weight"] = graph[:,3].tolist() # feel with free-flow travel times
    return g

#All_or_nothing assignment
def all_or_nothing(g, od):
    '''
    We are given an igraph object 'g' with od in the format {from: ([to], [rate])}
    do all_or_nothing assignment
    '''

    # csv to save ods that do not have paths in the graph
    csv_file = open("no_paths.csv", 'wb')
    writer = csv.writer(csv_file)
    writer.writerow(['orig', 'dest'])

    L = np.zeros(len(g.es), dtype="float64")
    count = 0

    for o in od.keys():

        out = g.get_shortest_paths(o, to=od[o][0], weights="weight", output="epath")

        for i, inds in enumerate(out):
            if len(inds) == 0:
                #print 'no path between {} and {}'.format(o, od[o][0][i])
                row = [o,od[o][0][i]]
                writer.writerow(row)
                count+=1
            L[inds] = L[inds] + od[o][1][i]
    csv_file.close()
    return L

def total_free_flow_cost(g, od):
    return np.array(g.es["weight"]).dot(all_or_nothing(g, od))

# Search directions step in Frank_Wolfe
def search_direction(f, bpr, g, od):
    # computes the Frank-Wolfe step
    # g is just a canvas containing the link information and to be updated with
    # the most recent edge costs
    x = np.power(f.reshape((f.shape[0],1)), np.array([0,1,2,3,4]))
    grad = np.einsum('ij,ij->i', x, bpr[:,3:])
    g.es["weight"] = grad.tolist()

    #start timer
    #start_time1 = timeit.default_timer()

    L = all_or_nothing(g, od)

    #end of timer
    #elapsed1 = timeit.default_timer() - start_time1
    #print ("all_or_nothing took  %s seconds" % elapsed1)

    return L, grad

#Calculating the potential of bpr function
def potential(graph ,f):
    # this routine is useful for doing a line search
    # computes the potential at flow assignment f
    links = int(np.max(graph[:,0])+1)
    g = graph.dot(np.diag([1.,1.,1.,1.,1/2.,1/3.,1/4.,1/5.]))
    x = np.power(f.reshape((links,1)), np.array([1,2,3,4,5]))
    return np.sum(np.einsum('ij,ij->i', x, g[:,3:]))

# Line Search step in Frank_Wolfe algorithm
def line_search(f, res=20):
    # on a grid of 2^res points bw 0 and 1, find global minimum
    # of continuous convex function
    d = 1./(2**res-1)
    l, r = 0, 2**res-1
    while r-l > 1:
        if f(l*d) <= f(l*d+d): return l*d
        if f(r*d-d) >= f(r*d): return r*d
        # otherwise f(l) > f(l+d) and f(r-d) < f(r)
        m1, m2 = (l+r)/2, 1+(l+r)/2
        if f(m1*d) < f(m2*d): r = m1
        if f(m1*d) > f(m2*d): l = m2
        if f(m1*d) == f(m2*d): return m1*d
    return l*d

def Frank_Wolfe_Solver(graph, demand, g=None, od=None, past=10, max_iter=100, eps=1e-16, \
    q=50, display=1, stop=1e-2):

    assert past <= q, "'q' must be bigger or equal to 'past'"
    if g is None:
        g = construct_igraph(graph)
    if od is None:
        od = construct_od(demand)
    f = np.zeros(graph.shape[0],dtype="float64") # initial flow assignment is null
    fs = np.zeros((graph.shape[0],past),dtype="float64") #not sure what fs does
    K = total_free_flow_cost(g, od)

    # why this?
    if K < eps:
        K = np.sum(demand[:,2])
    elif display >= 1:
        print ('average free-flow travel time', K / np.sum(demand[:,2]))

    for i in range(max_iter):

        if display >= 1:
            if i <= 1:
                print ('iteration: {}'.format(i+1))
            else:
                print ('iteration: {}, error: {}'.format(i+1, error))

        # construct weighted graph with latest flow assignment
        L, grad = search_direction(f, graph, g, od)

        fs[:,i%past] = L
        w = L - f
        if i >= 1:
            error = -grad.dot(w) / K
            # if error < stop and error > 0.0:
            if error < stop:
                if display >= 1: print ('stop with error: {}'.format(error))
                return f
        if i > q:
            # step 3 of Fukushima
            v = np.sum(fs,axis=1) / min(past,i+1) - f
            norm_v = np.linalg.norm(v,1)
            if norm_v < eps:
                if display >= 1: print ('stop with norm_v: {}'.format(norm_v))
                return f
            norm_w = np.linalg.norm(w,1)
            if norm_w < eps:
                if display >= 1: print ('stop with norm_w: {}'.format(norm_w))
                return f
            # step 4 of Fukushima
            gamma_1 = grad.dot(v) / norm_v
            gamma_2 = grad.dot(w) / norm_w
            if gamma_2 > -eps:
                if display >= 1: print ('stop with gamma_2: {}'.format(gamma_2))
                return f
            d = v if gamma_1 < gamma_2 else w
            # step 5 of Fukushima
            s = line_search(lambda a: potential(graph, f+a*d))
            if s < eps:
                if display >= 1: print ('stop with step_size: {}'.format(s))
                return f
            f = f + s*d
        else:
            f = f + 2. * w/(i+2.)

    return f

# Function to construct the od as dictionary of od and demand
def construct_od(demand):
    # construct a dictionary of the form
    # origin: ([destination],[demand])
    out = {}
    #import pdb; pdb.set_trace()
    for i in range(demand.shape[0]):
        origin = int(demand[i,0])
        if origin not in out.keys():
            out[origin] = ([],[])
        out[origin][0].append(int(demand[i,1]))
        out[origin][1].append(demand[i,2])
    return out

def main():
    # start timer for frank-wolfe
    start_time1 = timeit.default_timer()
    #Loading the graph data
    graph_data = np.loadtxt('bayarea_ter_igraph_bpr_coefficients.csv', delimiter=',', skiprows=1)
    #Loading the demand data
    demand = np.loadtxt('bayarea_ter_igraph_demand.csv', delimiter=',', skiprows=1)

    fileName = 'flow_on_Edges.csv'

    f = Frank_Wolfe_Solver(graph_data,demand)
    np.savetxt(fileName, f, delimiter=',')

    # end of timer
    elapsed1 = timeit.default_timer() - start_time1
    print ("Frank-Wolfe took  %s seconds" % elapsed1)



if __name__ == "__main__": main()