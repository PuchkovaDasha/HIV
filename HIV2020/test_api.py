import pytest

from api import Tree


def test_tree_creation():
    tree = Tree.from_vertex_distribution([2, 2, 2, 1])

    assert tree.levels[0].vertices == [0, 1]
    assert tree.levels[1].vertices == [2, 3]
    assert tree.levels[2].vertices == [4, 5]
    assert tree.levels[3].vertices == [6]


@pytest.mark.parametrize(
    "vertex_distribution",
    [
        [2, 2, 2, 1],
        [3, 3, 2, 1]
    ]
)
def test_random_correct_connections(vertex_distribution):
    tree = Tree.from_vertex_distribution(vertex_distribution)
    assert tree.random_correct_connections
    print(tree.random_correct_connections)
