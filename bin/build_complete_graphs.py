#!/usr/bin/env python3
import collections, os, csv, gzip, sys

# Disjoint-set data structure, without ranks
# Group 2-tuples if they share a parent node
# Partition of edges in connected components
class DisjointSet:
    def __init__(self):
        self.parent_pointer = collections.defaultdict(lambda: None)
        self.groups = collections.defaultdict(set)
    
    def find(self, x):
        curr_parent = self.parent_pointer[x]
        if curr_parent:
            updated_parent = self.find(curr_parent)
            self.parent_pointer[x] = updated_parent
            return updated_parent
        return x

    def union(self, x, y):
        parent_x, parent_y = self.find(x), self.find(y)
        if parent_x != parent_y:
            self.parent_pointer[parent_x] = parent_y
    
    def group(self):
        for x in self.parent_pointer:
            self.groups[self.find(x)].add(x)

# Fold: build k_n+1 from k_n
#   For each k_n
#    if exists v not in k_n: k_n U v = n and
#    for all u in k_n exists (u,v) then k_n+1 exists
def complete_graphs(vs: list, es: set) -> set:
    def build_kn(complete_graphs: set, n: int, acc: set) -> set:
      if n == len(vs) + 1: 
        return acc
      ks_step = set()
      for v in vs:
        for graph in complete_graphs:
            if len(graph | {v}) == n and all((v,u) in es or (u,v) in es for u in graph):
                ks_step |= {graph | {v}}
                acc |= {graph | {v}}
                acc -= {graph} - {frozenset({v})}
      return build_kn(ks_step, n+1, acc)
    return build_kn({frozenset([v]) for v in vs}, 1, set())

# Connect DisjointSet and complete graphs computations 
def complete_graphs_in_components(edges: list) -> list:
    disjoint_edges = DisjointSet()
    # Build tree
    for u,v in edges:
        disjoint_edges.union(str(u), str(v))
    # Build components
    disjoint_edges.group()
    components = list()
    for v in disjoint_edges.groups.values():
        v.discard(None)
        components.append(v)
    return [complete_graphs(component, edges) for component in components]

def cols_by_set(complete_graphs: list, prefixes: list) -> dict:
    acc, col = [], {prefix: [prefix] for prefix in prefixes}
    for component in complete_graphs:
        for graph in component:
            for k in col.keys():
                 l = [v for v in graph if v.startswith(k)]
                 assert len(l) <= 1, "Multiple vertices shared the same prefix"
                 if l:
                    col[k].append(l[0]) 
                 else:
                    col[k].append(None)
    return col

def cols_by_order(complete_graphs: list, prefixes: list) -> dict:
    n = len(prefixes)
    acc, col = [], {i: [i] for i in range(0,n)}
    for component in complete_graphs:
        for graph in component:
           l = [tuple(v.split("_")) for v in graph]
           # Definition of ordering, (+) i < j if i.start < j.start, (-) i < j if i.stop*-1 < j.stop*-1 
           l = sorted(l, key = lambda x: -1 * int(x[2]) if x[3] == "-" else int(x[1]))
           for i in range(0, n):
            try:
              col[i].append("_".join((str(s) for s in l[i])))
            except IndexError:
              col[i].append(None) 
    return col

# IO
ids, es, vs = set(), set(), set()
with gzip.open(sys.argv[1], mode='rt') as csvfile:
  rows = csv.reader(csvfile, delimiter=',', quotechar='|')
  for row in rows:
    i,j,max_distance,v,u,distance = row
    vs |= {v,u}
    ids |= {i,j}
    if abs(int(distance)) < abs(int(max_distance)):
      es |= {(v,u)}

# To return also k_1 graphs, loops are added
# Should slow down the computation and should be replaced by something more efficient
es.update([(vertice,vertice) for vertice in vs])
ids, gs = list(ids), complete_graphs_in_components(es)

for f in [cols_by_set, cols_by_order]:
   cols = f(gs, ids)
   with gzip.open(".".join([sys.argv[2], f.__name__, "csv", "gz"]), "wt") as out:
     wr = csv.writer(out)
     for x in zip(*[cols[k] for k in sorted(cols.keys())]): 
        wr.writerow(x) 
