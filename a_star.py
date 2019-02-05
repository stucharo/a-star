import heapq
from math import sin, cos, asin, pi, atan2, atan, acos
from copy import deepcopy
import numpy as np
from scipy.interpolate import griddata
from queue import PriorityQueue
import matplotlib.pyplot as plt
import random
from heuristic import shortest_route, copy_pt

spacing = 1
min_straight_length = 10
bend_radius = [2, 4, 6, 8, 10]
theta_tol = 3 * pi / 180
heading_tol = pi/180
location_tol = 1

class Grid:
    def __init__(self, eastings, northings, costs):
        self.eastings = eastings
        self.northings = northings
        self.costs = costs
        self.get_cost = self.create_interpolator()

    def create_interpolator(self):
        points = np.array((self.eastings.flatten(), self.northings.flatten())).T
        values = self.costs.flatten()
        def get_cost(X, Y):
            return griddata( points, values, (X,Y) )
        return get_cost

    def neighbours(self, pt):
        n = []
        if pt.radius != 0:
            # we keep on this bend or move the minimum straight length
            r = pt.radius
            a = spacing / r
            h = pt.heading + pt.direction * a
            p = copy_pt(pt, spacing*sin(h), spacing*cos(h), pt.direction*a)
            p.previous = pt
            p.cost = pt.cost + self.cost(pt, p)
            n.append(p)
            p = copy_pt(pt, min_straight_length*pt.dx, min_straight_length*pt.dy, 0)
            p.cost = pt.cost + self.cost(pt, p)
            n.append(p)
        else:
            # we can try every bend angle
            for r in bend_radius:
                a = spacing / r
                for i in [1, -1]:
                    h = pt.heading + i * a
                    new_x = pt.x+spacing*sin(h)
                    new_y = pt.y+spacing*cos(h)
                    p = Point(new_x, new_y, h, radius=)
                    p = copy_pt(pt, )
                    n.append()
                for f, d in [(1, 'r'), (-1, 'l')]:
                    h = pt.heading + f * a
                    new_x = pt.x+spacing*sin(h)
                    new_y = pt.y+spacing*cos(h)
                    p = Point(new_x, new_y, previous=pt, heading=h, radius=r)
                    p.cost = pt.cost + self.cost(pt, p)
                    n.append(p)
            new_x = pt.x+spacing*sin(pt.heading)
            new_y = pt.y+spacing*cos(pt.heading)
            p = copy_pt(pt, spacing*pt.dx, spacing*pt.dy, 0)
            p.previous = pt
            p.cost = pt.cost + self.cost(pt, p)
            n.append(p)
        return filter(self.is_valid, n)

    def is_valid(self, pt):
        return (np.amin(self.eastings) <= pt.x <= np.amax(self.eastings)
                and np.amin(self.northings) <= pt.y <= np.amax(self.northings))

    def cost(self, prev_p, p):
        d = dist(prev_p, p)
        pts = max(2, myround(d / spacing, spacing) + 1)
        xs = np.linspace(prev_p.x, p.x, pts)
        ys = np.linspace(prev_p.y, p.y, pts)
        costs = self.get_cost(xs, ys)
        cost = np.trapz(costs, dx=d/(pts-1))
        return cost
    
    def heuristic_cost(self, p, goal):
        path = shortest_route(p, goal, min(bend_radius), min_straight_length, heading_tol, location_tol, spacing)
        if len(path[0]) > 0:
            pts = np.array([(p.x, p.y) for p in path[0]])
            xs = pts[:,0]
            ys = pts[:,1]
            costs = self.get_cost(xs, ys)
            cost = np.trapz(costs, dx=spacing)
            return cost
        else:
            return 1_000_000

class Point:

    def __init__(self, x, y, heading, cost=0, radius=0, previous=None):
        self.x = x
        self.y = y
        self.dx = cos(heading)
        self.dy = sin(heading)
        self.cost = cost
        self.radius = radius
        self.previous = previous

    def __str__(self):
        return f"({self.x:.2f}, {self.y:.2f}) + ({self.dx:.2f}, {self.dy:.2f})"

    def __repr__(self):
        return f"({self.x:.2f}, {self.y:.2f}) + ({self.dx:.2f}, {self.dy:.2f})"

def plot_graph(start, goal, graph, path_ends):
    fig, ax = plt.subplots()
    ax.arrow(start.x, start.y, 5*sin(start.heading), 5*cos(start.heading), head_width=0.5, head_length=1)
    ax.arrow(goal.x, goal.y, 5*sin(goal.heading), 5*cos(goal.heading), head_width=0.5, head_length=1)
    for e in path_ends:
        p = e[2]
        pts = [(p.x, p.y)]
        while p.previous is not None:
            pts.append((p.previous.x, p.previous.y))
            p = p.previous
        pts = np.array(pts)
        plt.plot(pts[:,0], pts[:,1])
    plt.axis([0,100,0,100], option='equal')
    plt.show()

def a_star_pipe(graph, start, goal):
    path_ends = PriorityQueue()
    start.cost = 0
    path_ends.put((0, 0, start))
    counter = 1
    max_priority = 0

    while not path_ends.empty():
        # Get path end most likely to be cheapest
        end_point = path_ends.get()[2]

        if counter % 50 == 0:
            plot_graph(start, goal, graph, path_ends.queue)

        for neighbour in graph.neighbours(end_point):
            if neighbour == goal:
                return neighbour
            est_cost = neighbour.cost + graph.heuristic_cost(neighbour, goal)
            if np.isnan(est_cost) or np.isinf(est_cost):
                priority = 10 * max_priority
            else:
                priority = est_cost
                max_priority = max(max_priority, priority)
            print(f"{counter} {neighbour.x:.2f}, {neighbour.y:.2f}, {int(priority)}")
            path_ends.put((priority, counter, neighbour))
            counter += 1

    return end_point

if __name__ == '__main__':
    size = 100
    # create grid
    X = np.linspace(0, size, size+1)
    Y = np.linspace(0, size, size+1)
    X, Y = np.meshgrid(X, Y)
    Z = np.random.uniform(0,100, size=X.shape)
    graph = Grid(X, Y, Z)

    start = Point(10, 10, 0)
    goal = Point(90, 90, pi/2)

    pt = a_star_pipe(graph, start, goal)
    fig, ax = plt.subplots()
    ax.arrow(start.x, start.y, 5*sin(start.heading), 5*cos(start.heading), head_width=0.5, head_length=1)
    ax.arrow(goal.x, goal.y, 5*sin(goal.heading), 5*cos(goal.heading), head_width=0.5, head_length=1)
    pts = [(pt.x, pt.y)]
    while pt.previous is not None:
        pts.append((pt.previous.x, pt.previous.y))
        pt = pt.previous
    pts = np.array(pts)
    plt.plot(pts[:,0], pts[:,1])
    plt.axis([0,100,0,100], option='equal')
    plt.show()

