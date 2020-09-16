import itertools
import random
from dataclasses import dataclass
from itertools import groupby
from typing import List, Optional, Tuple

Vertex = int
Path = List[Vertex]
Connection = Tuple[Vertex, Vertex]


@dataclass
class Tree:
    vertices: List[Vertex]
    levels: List['Level']

    @classmethod
    def from_vertex_distribution(cls, vertex_distribution):
        vertices = list(range(sum(vertex_distribution)))
        levels = [
            Level([vertices.pop(0) for _ in range(level)])
            for level in vertex_distribution
        ]
        for level, next_level in zip(levels, levels[1:]):
            level.next_level = next_level

        vertices = tuple(range(sum(vertex_distribution)))
        return Tree(vertices, levels)

    @property
    def first_level_vertices(self):
        return self.levels[0].vertices

    @property
    def last_level_vertex(self):
        return self.levels[-1].vertices[0]

    @property
    def vertices_without_first_and_last_levels(self):
        return [
            vertex
            for vertex in self.vertices
            if vertex not in self.levels[0].vertices and
               vertex not in self.levels[-1].vertices
        ]

    @property
    def connections_generator(self):
        while True:
            connections = []

            for level in self.levels[:-1]:
                for vertex in level.vertices:
                    vertex_connections = level.list_vertex_connections(vertex)
                    connections.extend(
                        random.sample(vertex_connections, min(len(vertex_connections), 2))
                    )

            yield connections

    @property
    def duplication_part_vertices(self):
        return flat_map(
            lambda level: level.vertices,
            self.levels[1:1 + len(self.levels) // 3]
        )

    @property
    def random_correct_connections(self) -> List[Connection]:
        gen_ = self.connections_generator
        while True:
            connections = next(gen_)

            paths = Paths.from_connections(
                connections,
                self.first_level_vertices,
                self.vertices_without_first_and_last_levels
            ).paths

            all_paths_exist, path_cover_all_vertices = False, False
            for vertex_to_remove in self.duplication_part_vertices:
                paths_without_vertex = [path for path in paths if vertex_to_remove not in path.vertices]
                all_paths_exist = all(
                    any(vertex in path.vertices for path in paths_without_vertex)
                    for vertex in self.first_level_vertices
                )
                if not all_paths_exist:
                    break

                path_cover_all_vertices = all(
                    any(vertex in path.vertices for path in paths)
                    for vertex in self.vertices_without_first_and_last_levels
                )
                if not path_cover_all_vertices:
                    break

            if all_paths_exist and path_cover_all_vertices:
                return connections


@dataclass
class Paths:
    paths: List['Path']

    @classmethod
    def from_connections(
        cls, connections,
        first_level_vertices: List[Vertex],
        vertices_without_first_and_last_levels: List[Vertex]
    ) -> 'Paths':
        vertex_connections = {
            vertex: tuple(vertex_connections_)
            for vertex, vertex_connections_ in groupby(connections, lambda conn: conn[0])
        }

        path_connections = tuple(map(
            lambda connection: [connection],
            merge_lists(
                vertex_connections[vertex]
                for vertex in first_level_vertices
            )
        ))

        for vertex in vertices_without_first_and_last_levels:
            path_connections = flat_map(
                lambda connection: [
                    [*connection, c]
                    for c in vertex_connections[vertex]
                ] if vertex == connection[-1][-1] else
                [connection],
                path_connections
            )

        return Paths(paths=tuple(map(Path.from_connections, path_connections)))


@dataclass
class Path:
    vertices: List[Vertex]

    @classmethod
    def from_connections(cls, connections):
        """
        >>> Path.from_connections([[0,1], [1, 2], [2, 3]])
        Path(vertices=[0, 1, 2, 3])
        """
        return Path(sorted(frozenset(itertools.chain.from_iterable(connections))))


@dataclass
class Level:
    vertices: List[Vertex]
    next_level: Optional['Level'] = None

    def list_vertex_connections(self, vertex: Vertex):
        assert vertex in self.vertices
        return [
            (vertex, next_level_vertex)
            for next_level_vertex in self.next_level.vertices
        ]

    def connection_pairs(self):
        return tuple(
            itertools.chain.from_iterable(
                self.list_vertex_connections(vertex)
                for vertex in self.vertices
            )
        )


def merge_lists(list_of_lists):
    """
    >>> merge_lists([[1, 2], [3, 4]])
    [1, 2, 3, 4]
    """
    return tuple(itertools.chain.from_iterable(list_of_lists))


def flat_map(func, list_):
    """
    >>> flat_map(lambda item: [item, item], [1, 2])
    [1, 1, 2, 2]
    """
    return merge_lists(map(func, list_))
