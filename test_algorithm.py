from dataclasses import dataclass, field
from typing import Any, List
import random

import numpy as np

@dataclass
class Node:
    acc_cost: float
    p: Any = None
    open_n: List[int] = field(default_factory=list)
    best_n: Any = None
    r: float = None

    def __repr__(self):
        if self.best_n is None:
            bn = None
        else:
            bn = self.best_n.acc_cost
        return f"cost: {self.acc_cost:.2f}, r={self.r}, best_n: {bn}, open_ns: {len(self.open_n)}"

def get_route(s):
    cur = s
    for _ in range(random.randrange(1, 10)):
        n = Node(acc_cost=cur.acc_cost+random.uniform(100, 1000), r=random.choice([-5, 0, 5]))
        n.p = cur
        cur.best_n = n
        cur = n
    return cur

def print_route(g):
    cur = g
    while cur.p is not None:
        print(cur)
        cur = cur.p
    print(s)

def get_nexts(p):
    if p.r == 0:
        p.open_n.append(Node(acc_cost=p.acc_cost+random.uniform(100, 1000),r=-5, p=p))
        p.open_n.append(Node(acc_cost=p.acc_cost+random.uniform(100, 1000),r=5, p=p))
    else:
        p.open_n.append(Node(acc_cost=p.acc_cost+random.uniform(100, 1000),r=0, p=p))

def explore(s, best_cost, best_route):
    e = get_route(s)
    print(f"s.cost: {s.acc_cost:.2f}, e.cost: {e.acc_cost:.2f}")
    if e.acc_cost < best_cost:
        print_route(e)
        best_cost = e.acc_cost
        best_route = e
    cur = e.p
    while True:
        get_nexts(cur)
        while len(cur.open_n) > 0:
            n = cur.open_n.pop()
            if n.acc_cost < best_cost:
                explore(n, best_cost, best_route)
        # roll back a node
        if cur is s:
            break
        cur = cur.p
    return best_cost, best_route

# seed the pseudo-random number generator so our routine is always deterministic
random.seed(0)
best_cost = np.inf
best_route = None
# create our starting node
s = Node(acc_cost=0, r=0)

print(best_cost)
print(best_route)
# this is basically the start of our new algorithm
bc, br = explore(s, best_cost, best_route)
print(bc)
print(br)





