import random
import matplotlib.pyplot as plt
import numpy as np

## Replicating Tornberg 2018  
## The results for p_n != 0 do not match, so the relocation algorithm seems to be wrong...


class Node:
    """Class of Nodes"""
    def __init__(self):
        self.connection = []
        self.connection_num = []
        self.in_chamber = False
        self.activated = False

    def connect(self, node):
        self.connection.append(node)

    def disconnect(self, node):
        self.connection.remove(node)

    def activate(self):
        self.activated = True

    def list_edges(self, nodes):
        for connection in self.connection:
            self.connection_num.append(nodes.index(connection))
        return sorted(self.connection_num)


class Edge:
    """Class of Edges"""
    def __init__(self, node1, node2, nodes):
        self.lower = min(nodes.index(node1), nodes.index(node2))
        self.higher = max(nodes.index(node1), nodes.index(node2))


def create_edge(nodes, nodes_list, edges):
    """Create edge between nodes, try again if already exits"""
    while True:
        node1 = random.choice(nodes)
        node2 = random.choice(nodes)
        if (node1 != node2) and (node2 not in node1.connection):
            node1.connection.append(node2)
            node2.connection.append(node1)
            edges.append(Edge(node1, node2, nodes_list))
            break    


def create_network(N, k, c, p_n):
    """Create network with N nodes, c * N in-chamber nodes."""

    nodes = []
    edges = []

    # Create N nodes.
    for i in range(N):
        node = Node()
        nodes.append(node)

    # Create edges between nodes.
    for i in range(int(N * k / 2)):    
        create_edge(nodes, nodes, edges)

    # Determine the nodes inside echo chamber
    in_node =[]
    for node in random.sample(nodes, int(c*N)):
        node.in_chamber = True
        in_node.append(node)

    # Find the edges that span beyond clusters
    span_edges = []
    for edge in edges:
        if nodes[edge.higher].in_chamber != nodes[edge.lower].in_chamber:
            span_edges.append(edge)

    # Remove fraction of those edges
    relocate_num = int(p_n * len(span_edges))
    remove_list = random.sample(span_edges, relocate_num)
    for edge in remove_list:
        nodes[edge.lower].disconnect(nodes[edge.higher])
        nodes[edge.higher].disconnect(nodes[edge.lower])
        edges.remove(edge)

    # Relocate the same number within the echo chamber
    for i in range(relocate_num):        
        create_edge(in_node, nodes, edges)

    return [nodes, edges]


def activate_network(nodes, N, threshold, p_o):
    """Activate the nodes, return proportion activated"""

    in_node =[]

    # list of in-chamber nodes
    for node in nodes:
        if node.in_chamber:
            in_node.append(node)
    
    # Activate seed node and the surrounding nodes.
    seed = random.choice(in_node)
    seed.activate()
    for node in seed.connection:
        node.activate()

    # Count initially activated nodes.
    num_activated = 0
    for node in nodes:
        if node.activated:
            num_activated += 1
    num_activated_prev = 0

    # Activate the neighboring nodes until change stops
    while num_activated != num_activated_prev:
        num_activated_prev = num_activated

        for node in nodes:
            # ignore nodes with no connection
            if not node.activated and len(node.connection) != 0:
                activated_neighbor = 0
                for connection in node.connection:
                    if connection.activated:
                        activated_neighbor += 1
                frac = activated_neighbor / len(node.connection)
                # Differentiate threshold for in- and out-chamber nodes
                if (node.in_chamber and frac >= (threshold - p_o)) or (frac >= threshold):
                    node.activate()
                    num_activated += 1 
    
    return num_activated/N


def virality(N, k, c, p_n, threshold, p_o, test_num=1000):
    """Calculate the proportion of viral results"""

    viral = 0
    for i in range(test_num):
        network = create_network(N, k, c, p_n)[0]
        if activate_network(network, N, threshold, p_o) >= 0.5:
            viral += 1
    return viral/test_num


# Parameters
N = 100 # number of nodes
k = 8 # mean degree
c = 0.2 # proportion of in-chamber nodes
threshold = 0.27 # activation threshold
p_n = 0 # network polarization parameter
p_o = 0 # opinion polarization parameter

if __name__ == '__main__':
    """Run file, plot results"""
    p_n_list = np.arange(0, 1, 0.075)
    threshold_list = [0.150, 0.174, 0.198, 0.222, 0.246, 0.270, 
                        0.294, 0.318, 0.342, 0.366, 0.390, 0.414]
    
    for threshold in threshold_list:
        results = []
        for p_n in p_n_list:
            viral = virality(N, k, c, p_n, threshold, p_o)
            results.append(viral)
            print(viral)
        plt.plot(p_n_list, results)

    plt.xlabel('Network Polarization')
    plt.ylabel('Virality')
    plt.show()