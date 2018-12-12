import heapq
from math import sin, cos, asin, pi
from copy import deepcopy
import numpy as np
from scipy.interpolate import griddata
from queue import PriorityQueue
import matplotlib.pyplot as plt

from dubins import dubins_routes

spacing = 1
bend_radius = [5]
heading_tol = 3 * pi / 180

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
        angles = [spacing / r for r in bend_radius]
        angles.extend([-1 * a for a in angles])
        [n.append(pt.bend_pt(self, theta)) for theta in angles]
        # and single straight point
        n.append(pt.straight_pt(self, spacing))

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
        start = (p.x, p.y, p.heading)
        end = (goal.x, goal.y, goal.heading)
        r = dubins_routes(start, end, min(bend_radius), spacing)
        r.sort(key=lambda x: x[0])
        for p in r:
            xs = p[1][:,0]
            ys = p[1][:,1]
            costs = self.get_cost(xs, ys)
            cost = np.trapz(costs, dx=spacing)
            if not np.isnan(cost) and not np.isinf(cost):
                break
        return cost

def dist(p1, p2):
    return ((p1.x-p2.x)**2 + (p1.y-p2.y)**2)**0.5

def myround(num, div):
   whole = int(num/div)
   if (num%div)/div >= 0.5:
           return div * (whole + 1)
   else:
           return div * whole


class Point:
    """ Points are basically nodes in a linked list so we need to keep
    track of previous and subsequent Points. Each point can only have 1
    previous node (lying on the cheapest path to that point) but while
    the algorithm is exploring there may be several subsequent points so
    those should be stored in a list."""
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.theta = 0
        self.heading = 0
        self.cost = 0
        self.previous = None
    
    def bend_pt(self, graph, theta=None):
        if theta is None:
            new_heading = self.heading + self.theta
        else:
            new_heading = self.heading + theta
        new_x = self.x + spacing * sin(new_heading)
        new_y = self.y + spacing * cos(new_heading)
        p = Point(new_x, new_y)
        p.theta = self.theta
        p.heading = new_heading
        p.previous = self
        p.cost = self.cost + graph.cost(self, p)
        return p

    def straight_pt(self, graph, dist):
        new_x = self.x + dist * sin(self.heading)
        new_y = self.y + dist * cos(self.heading)
        p = Point(new_x, new_y)
        p.theta = 0
        p.heading = self.heading
        p.previous = self
        p.cost = self.cost + graph.cost(self, p)
        return p

    def __eq__(self, p):
        return dist(self, p) < 1*spacing and abs(self.heading - p.heading) < heading_tol

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
    plt.axis([0,20,0,20], option='equal')
    plt.show()

def a_star_pipe(graph, start, goal):
    path_ends = PriorityQueue()
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
            path_ends.put((priority, counter, neighbour))
            plt.plot([end_point.x, neighbour.x], [end_point.y, neighbour.y])
            counter += 1

    return end_point

# create grid
X = np.linspace(0, 20, 21)
Y = np.linspace(0, 20, 21)
X, Y = np.meshgrid(X, Y)
Z = np.random.uniform(0,100, size=X.shape)
graph = Grid(X, Y, Z)

start = Point(2, 2)
start.heading = 0
goal = Point(15, 15)
goal.heading = pi/2

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
plt.axis([0,20,0,20], option='equal')
plt.show()

