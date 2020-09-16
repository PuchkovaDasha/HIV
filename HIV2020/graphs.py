from abc import abstractmethod
from typing import List

import networkx as nx
from networkx.drawing.nx_agraph import to_agraph


class GraphvizGraphRender:
    def __init__(self, graph, filename):
        self.graph = graph
        self.filename = filename

    def render(self):
        # https://stackoverflow.com/questions/39657395/how-to-draw-properly-networkx-graphs
        A = to_agraph(self.graph)
        A.layout('dot')
        A.draw(self.filename)


class BaseGraphBuilder:
    colors = [
        "gold",
        "blue1",
        "green3"
    ]

    def __init__(self, tree):
        self.tree = tree

    @abstractmethod
    def create_nodes(self) -> List[dict]:
        ...

    @abstractmethod
    def create_edges(self) -> List[tuple]:
        ...

    def build(self) -> nx.DiGraph:
        G = nx.DiGraph()

        for node in self.create_nodes():
            G.add_node(**node)
        for edge in self.create_edges():
            G.add_edge(*edge)

        return G


class WithoutClosureGraphBuilder(BaseGraphBuilder):
    def create_nodes(self) -> List[dict]:
        nodes = []
        for index, level in enumerate(self.tree.levels):
            for vertex in level.vertices:
                nodes.append({
                    'node_for_adding': vertex,
                    'rank': index,
                    'color': self.colors[index % len(self.colors)]
                })
        return nodes

    def create_edges(self) -> List[tuple]:
        return self.tree.random_correct_connections


class WithClosureGraphBuilder(WithoutClosureGraphBuilder):
    @property
    def closure_node(self):
        return len(self.tree.vertices)

    def create_nodes(self) -> List[dict]:
        return [
            *super().create_nodes(),
            {'node_for_adding': self.closure_node, 'rank': len(self.tree.levels), 'color': 'red'}
        ]

    def create_edges(self) -> List[tuple]:
        return [
            *super().create_edges(),
            (self.closure_node - 1, self.closure_node),
            *[(self.closure_node, vertex) for vertex in self.tree.first_level_vertices]
        ]
