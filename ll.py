
from typing import Any, List
from dataclasses import dataclass, field

@dataclass
class RouteNodeChildList:
    """ A lightweight data class to act as a linked list of
    route nodes. Each node can only have one previous node (parent) but
    may have many route options from that point (child)
    """
    x: float
    y: float
    children: List[int]
    parent: int

    def __hash__(self):
        return hash(hash(self.x+self.y) + hash(self.parent))

points = {}

def children(p):
    return [points[c] for c in p.children]

def parent(p):
    return points[p.parent]


prev_nl = RouteNodeChildList(0, 0, [], None)

points[hash(prev_nl)] = prev_nl

for i in range(1, 10):
    nl = RouteNodeChildList(i, i, [], hash(prev_nl))
    points[hash(nl)] = nl
    prev_nl.children.append(hash(nl))
    prev_nl = nl

print(parent(nl))




