from queue import PriorityQueue
from math import pi
from dataclasses import dataclass
from typing import Any

import numpy as np
import matplotlib.pyplot as plt

from heuristic import RouteNode, bend_point, copy_pt, RouteNode, get_shortest_path, straight_points

@dataclass
class Route:
    cost: float
    loc: RouteNode = None

bend_rads = [5, 10]
min_straight = 10
heading_tol = pi/180
location_tol = 1
spacing = 1

"""
start at point
add heuristic cost to goal to current cost
if total is less than current best cost
    replace current best cost
move to parent node
add child to list of checked children
if other children to check:
    move to next point
    repeat
else:
    move to parent node
    repeat

"""

def cost(r):

    cost = 0
    pt = r
    while pt.prev is not None:
        cost += spacing
        pt = pt.prev

    return cost

def get_last_pt(p):
    cur = p
    while cur.child is not None:
        cur = cur.child
    return cur

def first_pt(p):

    while p.prev is not None:
        p = p.prev
    
    return p

def graph(path):

    ps = []
    while path.prev is not None:
        ps.append((path.x, path.y))
        path = path.prev
    ps.reverse()
    if len(ps) > 1:
        ps = np.array(ps)
        plt.plot(ps[:,0], ps[:,1])
        plt.axis('equal')
        plt.show()

def write_route(r):
    
    ps = []
    while r.prev is not None:
        ps.append((r.x, r.y))
        r = r.prev
    ps.reverse()
    for p in ps:
        print(p)

def alt_options(p):
    o = []
    s = p.prev
    if s.radius != 0:
        if p.radius != s.radius:
            o.append(bend_point(s, spacing, s.radius))
        else:
            g = RouteNode(s.x+min_straight*s.dx, s.y+min_straight*s.dy, s.dx, s.dy, 0)
            ns = straight_points(s, g, spacing)
            ns[0].prev = p
            for i in range(1, len(ns)):
                ns[i].prev = ns[i-1]
            o.append(ns[-1])
    else:
        d2b = 0
        _p = p
        while _p.prev.radius == 0:
            d2b += spacing
            _p = _p.prev
        if d2b > min_straight:
            neg_rads = [-r for r in bend_rads]
            rads = [r for r in bend_rads]
            rads.extend(neg_rads)
            o.append([bend_point(p, spacing, r) for r in rads])
    return o

def explore_backwards(s, g, best_route, best_cost):

    r = get_shortest_path(s, g, min(bend_rads),
        min_straight, min_straight, heading_tol, location_tol,
        spacing)
    c = cost(r)
    route = Route(c, r)
    
    if c < best_cost:
        best_cost = c
        best_route = route
    
    p = r

    while p.prev is not None:
        for opt in alt_options(p):
            if cost(opt) < best_cost:
                explore_backwards(opt, g, best_route, best_cost)
        p = p.prev
        

start = RouteNode(0, 0, 0, 1, 0, 0, None)
goal = RouteNode(100, 0, 0, -1, 0, 0, None)

best_route = None
best_cost = np.inf

explore_backwards(start, goal, best_route, best_cost)
     
