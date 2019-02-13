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
bend_radius = [10, 20]
theta_tol = 3 * pi / 180
heading_tol = pi/180
location_tol = 1

def neighbours(pt):
    n = []
    if pt.radius != 0:
        # we keep on this bend or move the minimum straight length
        n.append(bend_point(pt, spacing, pt.radius))
        sp = copy_pt(pt, min_straight_length*pt.dx, min_straight_length*pt.dy, 0)
        sp.radius = 0
        n.append(sp)
    else:
        # we can try every bend angle
        for r in bend_radius:
            for d in [1, -1]:
                dr = d * r
                n.append(bend_point(pt, spacing, dr))
        n.append(copy_pt(pt, spacing*pt.dx, spacing*pt.dy, 0))
    return n

def cost(grid, prev_p, p):
    d = dist(prev_p, p)
    pts = max(2, myround(d / spacing, spacing) + 1)
    xs = np.linspace(prev_p.x, p.x, pts)
    ys = np.linspace(prev_p.y, p.y, pts)
    costs = grid.get_cost(xs, ys)
    cost = np.trapz(costs, dx=d/(pts-1))
    return cost
    
def heuristic_cost(self, p, goal):
    pts = []
    cur_p = p
    while cur_p.previous is not None:
        pts.append((p.previous.x, p.previous.y))
        cur_p = cur_p.previous
    heuristic_path = get_shortest_path(p, goal, min(bend_radius), min_straight_length, min_straight_length, heading_tol, location_tol, spacing)
    #plot(pts + [(p.x, p.y) for p in heuristic_path])
    if len(heuristic_path) > 0:
        pts = np.array([(p.x, p.y) for p in heuristic_path])
        xs = pts[:,0]
        ys = pts[:,1]
        costs = self.get_cost(xs, ys)
        if np.isnan(costs.max()) or np.isnan(costs.min()):
            return 100_000_000
        cost = np.trapz(costs, dx=spacing)
        return cost
    else:
        return 100_000_000

def plot(path):
    p = np.asarray(path)
    xs = p[:,0]
    ys = p[:,1]
    plt.plot(xs, ys)
    plt.show()

def myround(num, div):
   whole = int(num/div)
   if (num%div)/div >= 0.5:
        return div * (whole + 1)
   else:
        return div * whole

def plot_graph(start, goal, graph, path_ends):
    fig, ax = plt.subplots()
    ax.arrow(start.x, start.y, 5*start.dx, 5*start.dy, head_width=0.5, head_length=1)
    ax.arrow(goal.x, goal.y, 5*goal.dx, 5*goal.dy, head_width=0.5, head_length=1)
    plt.contour(graph.eastings, graph.northings, graph.costs)
    for e in path_ends:
        p = e[2]
        pts = [(p.x, p.y)]
        while p.previous is not None:
            pts.append((p.previous.x, p.previous.y))
            p = p.previous
        pts = np.array(pts)
        plt.plot(pts[:,0], pts[:,1])
    plt.axis([0,50,0,50], option='equal')
    plt.show()

def a_star_pipe(start, goal):
    path_ends = PriorityQueue()
    start.cost = 0
    path_ends.put((0, 0, start))
    counter = 1
    max_priority = 0

    while not path_ends.empty():
        # Get path end most likely to be cheapest
        cur_point = path_ends.get()[2]

        neighbours = get_neighbours(cur_point)

        for neighbour in neighbours(end_point):
            neighbour.previous = end_point
            neighbour.cost = end_point.cost + cost(graph, end_point, neighbour)

            plot_graph(start, goal, graph, path_ends.queue)

            if neighbour == goal:
                return neighbour

            hp = get_shortest_path()
            hc = heuristic_cost(graph, neighbour, goal)

            if hc < np.inf:
                est_cost = neighbour.cost + hc
                print(f"{counter}: {100*neighbour.cost/est_cost:.2f}% route solved")
                if np.isnan(est_cost) or np.isinf(est_cost):
                    priority = 10 * max_priority
                else:
                    priority = est_cost
                    max_priority = max(max_priority, priority)
                path_ends.put((priority, counter, neighbour))
                counter += 1

    return end_point

if __name__ == '__main__':
    size = 50
    # create grid
    np.random.seed(0)
    X = np.linspace(0, size, size+1)
    Y = np.linspace(0, size, size+1)
    X, Y = np.meshgrid(X, Y)
    Z = np.ones(X.shape)
    graph = Grid(X, Y, Z)

    start = Point(10, 10, 0)
    goal = Point(40, 40, pi/2)

    pt = a_star_pipe(graph, start, goal)
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

