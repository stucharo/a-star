from queue import PriorityQueue
from math import pi
from dataclasses import dataclass
from typing import Any

import numpy as np

from heuristic import Point, bend_point, copy_pt, RouteNode, get_shortest_path

@dataclass
class Route:
    cost: float
    loc: RouteNode = None

bend_rads = [5, 10]
min_straight = 10
heading_tol = pi/180
location_tol = 1
spacing = 1

def cost(l):
    return l

def get_alt_routes(cr):
    # this returns all path options from the current
    # point EXCEPT the next point on the heuristic path
    ars = []
    l = cr.loc
    # if cp is on a bend
    if l.radius != 0:
        # we eith continue around the bend OR move straight
        if l.next.radius == 0:
            # return the bend point
            bp = bend_point(l, spacing, l.radius)
            bp.previous = cp
            bp.cost = bp.previous.cost + cost(spacing)
            aps.append(bp)
        else:
            #return the straight point the minimum straight distance away
            sp = copy_pt(cp, min_straight*cp.dx, min_straight*cp.dy, 0)
            sp.radius = 0
            sp.previous = cp
            sp.cost = sp.previous.cost + cost(min_straight)
            aps.append(sp)
    else:
        # we return all options except the heuristic option
        # first, lets make a list of all the bend options
        neg_rads = [-r for r in bend_rads]
        rads = [r for r in bend_rads]
        rads.extend(neg_rads)
        # lets check if the heuristic path moves into a bend
        if l.next.radius != 0:
            # remove that radius
            rads = list(filter(lambda x: x != l.next.radius), bend_rads)
            # add a straight point the a unit distance away
            new_l = copy_pt(l, spacing*l.dx, spacing*l.dy, 0)
            new_l.prev = l
            sr = Route(cr.cost + cost(spacing), new_l)
            ars.append(sr)
        # then add all the bend points
        for r in rads:
            bp = bend_point(l, spacing, r)
            bp.prev = l
            br = Route(cr.cost + cost(spacing), bp)
            ars.append(br)
    return ars

def count_pts(p):
    i = 1
    while p.next is not None:
        p = p.next
        i += 1
    return i

def route_finder(s, g):

    routes = PriorityQueue()
    counter = 0
    full_route = None
    # Estimate the first route
    r = Route(0, get_shortest_path(s, g, min(bend_rads), min_straight, min_straight, heading_tol, location_tol, spacing))
    # Move to the last tile (the goal)
    while r.loc.next is not None:
        r.loc = r.loc.next
    # calculate the total cost to 
    routes.put((cost(count_pts(r.loc.next) * spacing), counter, r))

    while not routes.empty():

        # get the priority path
        priority = routes.get()
        cr = priority[2]
        print(f"Current route: {cr.loc.x}, {cr.loc.y}")
        # move to the next point along the heuristic path
        cr.loc = cr.loc.next
        cr.cost += cost(spacing)
        if cr.loc.next.next is None:
            # the next point is the goal so there's no 
            # point in exploring any options off this path
            cr.loc = cr.loc.next
            cr.cost += cost(spacing)
            full_route = cr
            print(f"{cr.loc.x}, {cr.loc.y}")
            continue
        # and push that back onto the queue with the same priority
        counter += 1
        routes.put((cr.cost + cost(count_pts(cr.loc.next)), counter, cr))

        # now get alternative path options from this point
        alt_routes = get_alt_routes(cr)
        for r in alt_routes:
            r.next = get_shortest_path(r.loc, g, min(bend_rads), min_straight, min_straight, heading_tol, location_tol, spacing).next
            if count_pts(r.next) > 0:
                priority = r.cost+cost(count_pts(r.next)*spacing)
                if full_route is not None and r.cost < full_route.cost:
                    counter += 1
                    routes.put((priority, counter, r))
    
    return full_route

start = RouteNode(10, 10, 1, 0, 0)
goal = RouteNode(40, 40, 0, 1, 0)

pt = route_finder(start, goal)
print(pt)


        
