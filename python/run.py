import networkx as nx
import plotly.graph_objs as go
import matplotlib.pyplot as plt

arcs = [ 
 5 -> 33;
 30 -> 7;
 30 -> 2;
 32 -> 45;
 32 -> 41;
 32 -> 19;
 32 -> 8;
 6 -> 33;
 45 -> 41;
 45 -> 19;
 45 -> 8;
 4 -> 33;
 13 -> 52;
 13 -> 47;
 13 -> 14;
 12 -> 13;
 12 -> 52;
 12 -> 28;
 12 -> 11;
 12 -> 26;
 12 -> 27;
 12 -> 47;
 12 -> 14;
 12 -> 51;
 12 -> 33;
 28 -> 13;
 28 -> 52;
 28 -> 11;
 28 -> 26;
 28 -> 27;
 28 -> 47;
 28 -> 14;
 23 -> 30;
 23 -> 7;
 23 -> 50;
 23 -> 2;
 23 -> 42;
 23 -> 20;
 41 -> 8;
 7 -> 2;
 25 -> 13;
 25 -> 52;
 25 -> 12;
 25 -> 28;
 25 -> 11;
 25 -> 26;
 25 -> 27;
 25 -> 47;
 25 -> 14;
 25 -> 51;
 25 -> 33;
 50 -> 30;
 50 -> 7;
 50 -> 2;
 50 -> 20;
 2 -> 7;
 10 -> 43;
 10 -> 9;
 18 -> 7;
 18 -> 2;
 18 -> 17;
 26 -> 13;
 26 -> 52;
 26 -> 28;
 26 -> 11;
 26 -> 27;
 26 -> 47;
 26 -> 14;
 27 -> 13;
 27 -> 52;
 27 -> 28;
 27 -> 11;
 27 -> 26;
 27 -> 47;
 27 -> 14;
 42 -> 7;
 42 -> 2;
 20 -> 30;
 20 -> 7;
 20 -> 50;
 20 -> 2;
 19 -> 41;
 19 -> 8;
 31 -> 7;
 31 -> 2;
 31 -> 42;
 31 -> 21;
 31 -> 17;
 21 -> 7;
 21 -> 2;
 21 -> 42;
 21 -> 17;
 17 -> 7;
 17 -> 2;
 37 -> 38;
 37 -> 48;
 47 -> 13;
 47 -> 52;
 47 -> 14;
 14 -> 52;
 3 -> 41;
 3 -> 19;
 3 -> 8;
 51 -> 11;
 51 -> 33;
 40 -> 38;
 40 -> 48;
 48 -> 38;
 15 -> 33;
 ]

D = nx.DiGraph()

for a in arcs:
    D.add_edge(a[0], a[1])

pos = nx.spring_layout(D, k=0.5, iterations=100)
for n, p in pos.items():
    D.nodes[n]['pos'] = p

# arcs
edge_trace = go.Scatter(
        x=[],
        y=[],
        marker=dict(size=15, symbol="arrow-bar-up", angleref="previous"),
        )

for edge in D.edges():
    x0, y0 = D.nodes[edge[0]]['pos']
    x1, y1 = D.nodes[edge[1]]['pos']
    edge_trace['x'] += tuple([x0, x1, None])
    edge_trace['y'] += tuple([y0, y1, None])

# nodes
node_trace = go.Scatter(
        x=[],
        y=[],
        text=[],
        mode='markers+text',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            colorscale='pinkyl',
            reversescale=True,
            color=[],
            size=37,
            colorbar=dict(
                thickness=1,
                title='Node Connections',
                xanchor='left',
                titleside='right'
                ),
            line=dict(width=0)
            )
        )

for node in D.nodes():
    x, y = D.nodes[node]['pos']
    node_trace['x'] += tuple([x])
    node_trace['y'] += tuple([y])

    for node, adjacencies in enumerate(D.adjacency()):
        node_trace['marker']['color'] += tuple([len(adjacencies[1])])
        node_info = adjacencies[0]
        node_trace['text'] += tuple([node_info])

#
title = "Network Graph Demonstration"
fig = go.Figure(data=[edge_trace, node_trace],
                layout=go.Layout(
                    title=title,
                    titlefont=dict(size=16),
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=21, l=5, r=5, t=40),
                    annotations=[dict(
                        text="Text Here",
                        showarrow=False,
                        xref="paper", yref="paper")],
                    xaxis=dict(showgrid=False, zeroline=False,
                               showticklabels=False, mirror=True),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, mirror=True)))
fig.show()

#nx.draw_networkx_nodes(D, pos)
#nx.draw_networkx_labels(D, pos)
#nx.draw_networkx_edges(D, pos, edge_color='r', arrows = True)
#plt.show()

"""
G = nx.Graph()

for a in arcs:
    G.add_edges_from([a])

pos = nx.spring_layout(G, k=0.5, iterations=100)
for n, p in pos.items():
    G.nodes[n]['pos'] = p

edge_trace = go.Scatter(
        x=[],
        y=[],
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

for edge in G.edges():
    x0, y0 = G.nodes[edge[0]]['pos']
    x1, y1 = G.nodes[edge[1]]['pos']
    edge_trace['x'] += tuple([x0, x1, None])
    edge_trace['y'] += tuple([y0, y1, None])

node_trace = go.Scatter(
        x=[],
        y=[],
        text=[],
        mode='markers+text',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            colorscale='pinkyl',
            reversescale=True,
            color=[],
            size=37,
            colorbar=dict(
                thickness=1,
                title='Node Connections',
                xanchor='left',
                titleside='right'
                ),
            line=dict(width=0)
            )
        )
for node in G.nodes():
    x, y = G.nodes[node]['pos']
    node_trace['x'] += tuple([x])
    node_trace['y'] += tuple([y])

    for node, adjacencies in enumerate(G.adjacency()):
        node_trace['marker']['color'] += tuple([len(adjacencies[1])])
        node_info = adjacencies[0]
        node_trace['text'] += tuple([node_info])

title = "Network Graph Demonstration"
fig = go.Figure(data=[edge_trace, node_trace],
                layout=go.Layout(
                    title=title,
                    titlefont=dict(size=16),
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=21, l=5, r=5, t=40),
                    annotations=[dict(
                        text="Text Here",
                        showarrow=False,
                        xref="paper", yref="paper")],
                    xaxis=dict(showgrid=False, zeroline=False,
                               showticklabels=False, mirror=True),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, mirror=True)))
fig.show()
"""
