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

def graph(r):

    p = r
    pts = []
    while p.prev is not None:
        pts.append((p.x, p.y))
        p = p.prev
    
    if len(pts) > 0:
        pts.reverse()
        pts = np.array(pts)
        plt.plot(pts[:,0], pts[:,1])
        plt.axis('equal')
        plt.show()




def minimise_cost(s, g, best_route, best_cost):
    r = get_shortest_path(s, g, min(bend_rads),
        min_straight, min_straight, heading_tol, location_tol,
        spacing)
    if r is not None:
        graph(r)
        add_costs(r)
        print(s.x, s.y, g.x, g.y)
        if r.cost < best_cost:
            best_route = r
            best_cost = r.cost
        else:
            while r.prev.cost >= best_cost and r.prev is not s:
                r = r.prev
        # here we're sitting at the node AHEAD of the node we want to explore from
        # whether it's the cheapest route or we've rewound until we're low enough
        # cost again. Either way, we want to explore alternatives to the the previous
        # node.
        while r.prev not in [None, s]:
            r = r.prev
            alt_opts = alt_options(r)
            while len(alt_opts) > 0:
                minimise_cost(alt_opts.pop(0), g, best_route, best_cost)

def add_costs(r):

    cost = 0
    pt = r
    # rewind to start
    while pt.prev is not None and pt.cost == 0:
        last_pt = pt
        pt = last_pt.prev
    pt.next = last_pt
    while pt.next is not None:
        pt.next.cost = pt.cost + spacing
        pt = pt.next

def alt_options(p):
    # here we get alternative nodes to THIS NODE. That means, go back to the node before
    # this one and list the other options
    prev = p.prev
    o = []
    if prev.radius != 0:
        # the previous node was on a bend
        if p.radius == 0:
            # the option we're finding alternatives to WASN'T the continuation of a bend so
            # we can add the bend alternative
            o.append(bend_point(prev, spacing, prev.radius))
        else:
            # the option we're finding alternatives to WAS the continuation of a bend so
            # we can add the minimum straight alternative
            g = RouteNode(prev.x+min_straight*prev.dx, prev.y+min_straight*prev.dy, prev.dx, prev.dy, 0)
            ns = straight_points(prev, g, spacing)[0]
            ns[0].prev = prev  #check that the first straight point is the next point after our current point
            prev.next = ns[0]
            for i in range(1, len(ns)):
                ns[i].prev = ns[i-1]
                ns[i-1].next = ns[i]
            o.append(ns[-1])
    else:
        # we're on a straight section....are we far enough along to let us start a new bend
        d2b = 0
        _p = prev
        # count the spacing back to the previous bend
        while _p.prev.radius == 0:
            d2b += spacing
            _p = _p.prev
        if d2b > min_straight:
            neg_rads = [-r for r in bend_rads]
            rads = [r for r in bend_rads]
            rads.extend(neg_rads)
            o.extend([bend_point(prev, spacing, r) for r in rads])
    return o


start = RouteNode(0, 0, 0, 1, 0, 0, None)
goal = RouteNode(100, 0, 0, -1, 0, 0, None)

best_route = None
best_cost = np.inf

minimise_cost(start, goal, best_route, best_cost)

