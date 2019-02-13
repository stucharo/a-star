import heapq
from math import sin, cos, asin, pi, atan2, atan, acos
from copy import deepcopy
import numpy as np
from scipy.interpolate import griddata
from queue import PriorityQueue
import matplotlib.pyplot as plt
import random
from heuristic import Point, copy_pt, dist, get_shortest_path, bend_point, direction

spacing = 1
min_straight_length = 10
bend_radius = [5]
theta_tol = 3 * pi / 180
heading_tol = pi/180
location_tol = 1

def get_neighbours(pt, goal):
    ns = []
    if pt.radius != 0:
        # if we are on a bend and approaching the intersection with the last straight, we need to make sure that the bend neighbour doesn't continue further around the bend than neccessary and leave us unable to get to the end vector
        # we keep on this bend or move the minimum straight length
        bp = bend_point(pt, spacing, pt.radius)
        bp.previous = pt
        bp.cost = bp.previous.cost + cost(spacing)
        ns.append(bp)
        sp = copy_pt(pt, min_straight_length*pt.dx, min_straight_length*pt.dy, 0)
        sp.radius = 0
        sp.previous = pt
        sp.cost = sp.previous.cost + cost(min_straight_length)
        ns.append(sp)
    else:
        # we can try every bend angle
        for r in bend_radius:
            for d in [1, -1]:
                dr = d * r
                bp = bend_point(pt, spacing, dr)
                bp.previous = pt
                bp.cost = bp.previous.cost + cost(spacing)
                ns.append(bp)
        sp = copy_pt(pt, spacing*pt.dx, spacing*pt.dy, 0)
        sp.previous = pt
        sp.cost = sp.previous.cost + cost(spacing)
        ns.append(sp)
    return ns

def cost(l):
    return l

def a_star_pipe(start, goal):
    path_ends = PriorityQueue()
    start.cost = 0
    path_ends.put((0, 0, start))
    counter = 1
    max_priority = 0

    while not path_ends.empty():
        # Get path end most likely to be cheapest
        cur_point = path_ends.get()[2]

        neighbours = get_neighbours(cur_point, goal)

        for neighbour in neighbours:
            
            if dist(neighbour, goal) < location_tol and abs(neighbour.heading - goal.heading) < heading_tol:
                return neighbour

            hp = neighbour.get_heuristic_path(goal, min(bend_radius),
                min_straight_length, min_straight_length, heading_tol,
                location_tol, spacing)
            
            #graph(neighbour, goal, hp)
            if len(hp) == 0:
                hc = 100_000_000
            else:
                graph(neighbour, goal, hp)
                hc = cost(len(hp) * spacing)

            priority = neighbour.cost + hc

            print(f"{counter}: {100*neighbour.cost/priority:.2f}% route solved")

            path_ends.put((priority, counter, neighbour))

            counter += 1

    return end_point

def graph(s, g, path):

    c = s
    prev_path = []
    while c.previous is not None:
        prev_path.append(c.previous)
        c = c.previous
    prev_path.reverse()
    prev_path.extend(path)
    if len(path) > 1:
        ps = np.array([(p.x, p.y) for p in prev_path])
        plt.plot(ps[:,0], ps[:,1])
        plt.axis('equal')
        plt.title(f"Cost: {cost(len(ps) * spacing)}")
    else:
        cx = (s.x+g.x)/2
        cy = (s.y+g.y)/2
        w = max(abs(g.x-s.x), abs(g.y-s.y))
        plt.axis([cx-w/2-10,cx+w/2+10,cy-w/2-10,cy+w/2+10])
    plt.arrow(s.x, s.y, 5*s.dx, 5*s.dy, head_width=0.5, color='green')
    plt.arrow(g.x, g.y, 5*g.dx, 5*g.dy, head_width=0.5, color='red')
    plt.show()

if __name__ == '__main__':

    start = Point(10, 10, 0)
    goal = Point(40, 40, pi/2)

    pt = a_star_pipe(start, goal)
    fig, ax = plt.subplots()
    ax.arrow(start.x, start.y, 5*start.dx, 5*start.dy, head_width=0.5, head_length=1)
    ax.arrow(goal.x, goal.y, 5*goal.dx, 5*goal.dy, head_width=0.5, head_length=1)
    pts = [(pt.x, pt.y)]
    while pt.previous is not None:
        pts.append((pt.previous.x, pt.previous.y))
        pt = pt.previous
    pts = np.array(pts)
    plt.plot(pts[:,0], pts[:,1])
    plt.axis([0,100,0,100], option='equal')
    plt.show()

