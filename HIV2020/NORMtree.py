import random
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from networkx.drawing.nx_agraph import to_agraph

from api import Tree

# Пстроение графа
from HIV2020.graphs import WithoutClosureGraphBuilder, GraphvizGraphRender, WithClosureGraphBuilder

print("Введите кол-во узлов")
N = int(input())  # кол-во узлов
# N = 9

names_of_nodes = list(range(1, N + 1))  # названия узлов
print("названия узлов")
print(names_of_nodes)

# колво узлов первого уровня
N_1 = N // 3

print('колво узлов первого уровня  = ', N_1)

# 3. Определить колво узлов в остальных уровнях кроме первого и терминального
N_levels = N - N_1 - 1  # колво узлов в остальных уровнях кроме первого и терминального
print("колво узлов в остальных уровнях кроме первого и терминального", N_levels)

# 4. Определить число уровней
delit = N_levels // N_1
if N_levels % N_1 == 0:
    # мин число уровней, если упакованы плотно, мах
    k = random.randint(delit, N_levels // 2)
else:
    k = random.randint(delit + 1, N_levels // 2)
print('число среднх уровней равно k= ', k)

# 5. Определить колво узлов в каждом уровне:

#INPUT_NUMBER = 13
#MAX_PART = 6
#PARTS_CNT = 4

parts = [1, ] * k
result_parts = []
for _ in range(N_levels - k):
    random_part = random.randint(0, len(parts) - 1)
    parts[random_part] += 1
    if parts[random_part] >= N_1:
        result_parts.append(parts.pop(random_part))

diff_dop = sorted(result_parts + parts, reverse=True)

print("колво узлов в средних узлах diff_dop", diff_dop)

diff_dop.append(1)
diff_dop.insert(0, N_1)
print("число узлов во всех рядах diff_dop", diff_dop)

print()

all_levels = 2 + k
print('общее колво уровней', all_levels)

tree = Tree.from_vertex_distribution(vertex_distribution=diff_dop)

builder = WithoutClosureGraphBuilder(tree)
closure_builder = WithClosureGraphBuilder(tree)

GraphvizGraphRender(builder.build(), "graph_without_closure.png").render()
GraphvizGraphRender(closure_builder.build(), "graph_with_closure.png").render()

connections = closure_builder.create_edges()

print(connections) # cвязи

W = np.zeros([N + 1,N +1])
for i in range(N + 1):
    for j in range(N + 1):
        if (i, j) in connections:
            W[i][j] = 1

print(W)

#посчитать степени вершин в графе
stepen = [0 for i in range(N+1)]
for i in range (N+1):
    for j in range (N+1):
        stepen[i] += W[i][j]
        stepen[i] += W[j][i]

print('cтепень',stepen)

f = open('stepen_vershin.txt','a')
f.write(", ".join(map(str, stepen))+'\n')
f.close()

# pos = nx.spring_layout(G)
# nx.draw(G,pos, node_color = [*"g" * len(tree.vertices), "r"], node_size=800,with_labels=True)
# plt.show()