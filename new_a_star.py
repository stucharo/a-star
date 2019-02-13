from queue import PriorityQueue
from math import pi

import numpy as np

from heuristic import Point, bend_point

bend_rads = [5, 10]
min_straight = 10
heading_tol = pi/180
location_tol = 1
spacing = 1

def cost(l):
    return l

def get_alt_points(cp):
    # this returns all path options from the current
    # point EXCEPT the next point on the heuristic path
    aps = []
    # if cp is on a bend
    if cp.radius != 0:
        # we eith continue around the bend OR move straight
        if cp.heuristic_path[0].radius == 0:
            # return the bend point
            bp = bend_point(cp, spacing, cp.radius)
            bp.previous = cp
            bp.cost = bp.previous.cost + cost(spacing)
            aps.append(bp)
        else:
            #return the straight point the minimum straight distance away
            sp = copy_pt(cp, min_straight*pt.dx, min_straight*pt.dy, 0)
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
        if cp.heuristic_path[0].radius != 0:
            # remove that radius
            rads = list(filter(lambda x: x != cp.heuristic_path[0].radius), bend_rads)
            # add a straight point the a unit distance away
            sp = copy_pt(cp, spacing*cp.dx, spacing*cp.dy, 0)
            sp.radius = 0
            sp.previous = cp
            sp.cost = sp.previous.cost + cost(spacing)
            aps.append(sp)
        # then add all the bend points
        for r in rads:
            bp = bend_point(cp, spacing, r)
            bp.previous = cp
            bp.cost = bp.previous.cost + cost(spacing)
            aps.append(bp)
    return aps

def a_star(s, g):

    paths = PriorityQueue()
    counter = 0
    route_cost = np.inf
    last_tile = None
    s.get_heuristic_path(g, min(bend_rads), min_straight, min_straight, heading_tol, location_tol, spacing)
    paths.put((s.cost + cost(len(s.heuristic_path) * spacing), counter, s))

    while not paths.empty():

        # get the priority path
        priority = paths.get()
        current_point = priority[2]
        print(f"Current point: {current_point}")
        # move to the next point along the heuristic path
        next_point = current_point.heuristic_path[0]
        next_point.previous = current_point
        next_point.cost = current_point.cost + cost(spacing)
        if len(current_point.heuristic_path) == 1:
            # the next point is the goal so there's no 
            # point in exploring any options off this path
            route_cost = next_point.cost
            last_tile = next_point
            print(last_tile)
            continue
        else:
            next_point.heuristic_path = current_point.heuristic_path[1:]
        # and push that back onto the queue with the same priority
        counter += 1
        paths.put((priority[0], counter, next_point))

        # now get alternative path options from this point
        alt_points = get_alt_points(current_point)
        for p in alt_points:
            p.get_heuristic_path(g, min(bend_rads), min_straight, min_straight, heading_tol, location_tol, spacing)
            if len(p.heuristic_path) > 0:
                priority = p.cost+cost(len(p.heuristic_path)*spacing)
                if priority < route_cost:
                    counter += 1
                    paths.put((p.cost+cost(len(p.heuristic_path)*spacing),
                               counter, p))
    
    return route_cost, last_tile

start = Point(10, 10, 0)
goal = Point(40, 40, pi/2)

pt = a_star(start, goal)
print(pt)


        
