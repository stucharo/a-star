from dataclasses import dataclass
from typing import Any

@dataclass
class Node:
    x: int
    y: int
    next: Any = None
    prev: Any = None

def list_to_cls():

    pts = [(0,0), (1,1), (2,2), (3,3)]

    s = Node(pts[0][0], pts[0][1])
    prev = s

    for p in pts[1:]:
        n = Node(p[0], p[1], None, prev)
        prev.next = n
        prev = n

    return s