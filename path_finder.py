from queue import PriorityQueue
from math import pi
from dataclasses import dataclass
from typing import Any

import numpy as np
import matplotlib.pyplot as plt

from heuristic import RouteNode, bend_point, copy_pt, RouteNode, get_shortest_path, straight_points

bend_rads = [5, 10]
min_straight = 10
heading_tol = pi/180
location_tol = 1
spacing = 1


def minimise_cost(s, g, best_route, best_cost):
    print(s.x, s.y)
    r = get_shortest_path(s, g, min(bend_rads),
        min_straight, min_straight, heading_tol, location_tol,
        spacing)
    if r is not None:
        add_costs(r)
        if r.cost < best_cost:
            best_route = r
            best_cost = r.cost
        else:
            while r.prev.cost > best_cost:
                r = r.prev
        while r.prev not in [None, s]:
            r = r.prev
            alt_opts = alt_options(r)
            while len(alt_opts) > 0:
                minimise_cost(alt_opts.pop(0), g, best_route, best_cost)

def add_costs(r):

    cost = 0
    pt = r
    # rewind to start
    while pt.prev is not None:
        pt = pt.prev
    while pt.next is not None:
        pt.next.cost = pt.cost + spacing
        pt = pt.next

def alt_options(p):
    o = []
    s = p.prev
    if s.radius != 0:
        if p.radius != s.radius:
            o.append(bend_point(s, spacing, s.radius))
        else:
            g = RouteNode(s.x+min_straight*s.dx, s.y+min_straight*s.dy, s.dx, s.dy, 0)
            ns = straight_points(s, g, spacing)[0]
            ns[0].prev = s
            s.next = ns[0]
            for i in range(1, len(ns)):
                ns[i].prev = ns[i-1]
                ns[i-1].next = ns[i]
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
            o.extend([bend_point(p, spacing, r) for r in rads])
    if len(o) > 0 and type(o[0]) is list:
        breakpoint()
    return o


start = RouteNode(0, 0, 0, 1, 0, 0, None)
goal = RouteNode(100, 0, 0, -1, 0, 0, None)

best_route = None
best_cost = np.inf

minimise_cost(start, goal, best_route, best_cost)

