import networkx as nx
import graphviz
from networkx.drawing.nx_agraph import graphviz_layout
import copy
import pandas as pd
import numpy as np
import ast

def compute_graph(pens_addr):
    pens = pd.read_csv(pens_addr)

    axes = ["cohort1", 
            "cohort2", 
            "cohort3", 
            "cohort4"]
    
    options = [[[True], [False]],
               [[True], [False]],
               [[True], [False]],
               [[True], [False]]]
    
    red = (0x9C, 0x11, 0x0C)
    blue = (0x0C, 0x0E, 0x9C)
    white = (0xFF, 0xFF, 0xFF)
    def pen2color(pen, hi, lo):
        r_val = red[0] * (hi - pen) / (hi - lo) + white[0] * (pen - lo) / (hi - lo)
        g_val = red[1] * (hi - pen) / (hi - lo) + white[1] * (pen - lo) / (hi - lo)
        b_val = red[2] * (hi - pen) / (hi - lo) + white[2] * (pen - lo) / (hi - lo)
        return "#" + "{:02x}".format(int(r_val)) \
                   + "{:02x}".format(int(g_val)) \
                   + "{:02x}".format(int(b_val))
    
    start = [[True, False],
             [True, False],
             [True, False],
             [True, False]]
    
    layers = [[start]]
    edges = []
    
    for i in range(len(axes)):
        new_layer = []
        for j in range(len(layers[-1])):
            for k in range(len(options[i])):
                node = copy.deepcopy(layers[-1][j])
                node[i] = options[i][k]
                new_layer.append(node)
                edges.append((str(layers[-1][j]), str(node)))
        layers.append(new_layer)
        
    def node2filter(node):
      c1, c2, c3, c4 = node
      filtered = pens[pens["cohort1"].isin(c1) & \
                      pens["cohort2"].isin(c2) & \
                      pens["cohort3"].isin(c3) & \
                      pens["cohort4"].isin(c4) & \
                      pens["mode"].isin(["map"])]
      
      return filtered
    
    def filter2pen(filtered):
        return np.mean(filtered["penetrance"])
        
    possible_pens = []
    for layer in layers:
        for node in layer:
            possible_pens.append(filter2pen(node2filter(node)))
    hi = max(possible_pens)
    lo = min(possible_pens)
    
    tree = nx.DiGraph()
    for layer in layers:
        for node in layer:
            label = str(node).replace("], [", "]\n [")[1:-1]
            filtered = node2filter(node)
            if filtered.size == 0:
                tree.add_nodes_from([(str(node), {"label": label, "color": "#888888"})])
            else:
                color = pen2color(filter2pen(node2filter(node)), hi, lo)
                tree.add_nodes_from([(str(node), {"label": label, "color": color})])
            
    for edge in edges:
        tree.add_edge(edge[0], edge[1])
    
    positions = graphviz_layout(tree, prog='dot')
    
    nodes = []
    node_x = []
    node_y = []
    colors = []
    for node in tree.nodes:
        nodes.append(node)
        x, y = positions[node]
        node_x.append(x)
        node_y.append(y)
        colors.append(tree.nodes[node]["color"])
    
    edge_x = []
    edge_y = []
    for edge in tree.edges():
        x0, y0 = positions[edge[0]]
        x1, y1 = positions[edge[1]]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)
        
    return {"nodes": nodes,
            "penetrances": possible_pens,
            "node_x": node_x, 
            "node_y": node_y, 
            "colors": colors, 
            "edge_x": edge_x, 
            "edge_y": edge_y}
        
compute_penetrance_tree(node, lvwt_th, sbp_th, quadrant, pens_addr):
    c1, c2, c3, c4 = ast.literal_eval(node)
    pens = pd.read_csv(pens_addr)

    filtered = pens[
      pens["cohort1".isin(c1)] & \
      pens["cohort2".isin(c2)] & \
      pens["cohort3".isin(c3)] & \
      pens["cohort4".isin(c4)] & \
      pens["mode"].isin(["map"])
    ]
    
    if lvwt_th > 0:
        filtered = filtered[
          filtered["LVWT threshold"] == lvwt_th
        ]
    if sbp_th > 0:
        filtered = filtered[
          filtered["SBP threshold"] == lvwt_th
        ]
    if quadrant > 0:
        filtered = filtered[
          filtered["Quadrant"] == quadrant
        ]
    return np.mean(filtered["penetrance"])
