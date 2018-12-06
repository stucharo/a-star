from collections import namedtuple
import heapq
from math import sin, cos, asin, pi
from copy import deepcopy
import numpy as np
from scipy.interpolate import griddata

class PriorityQueue:
    def __init__(self):
        self.elements = []
    
    def empty(self):
        return len(self.elements) == 0
    
    def put(self, item, priority):
        heapq.heappush(self.elements, (priority, item))
    
    def get(self):
        return heapq.heappop(self.elements)[1]

def reconstruct_path(came_from, start, goal):
    current = goal
    path = []
    while current != start:
        path.append(current)
        current = came_from[current]
    path.append(start) # optional
    path.reverse() # optional
    return path

def heuristic(a, b):
    return ((a.x-b.x)**2 + (a.y-b.y)**2)**0.5

spacing = 1
min_straight_length = 10
min_bend_radius = 10
max_bend_radius = 200
bend_radius_inc = 10
heading_tol = pi / 180

class Grid:
    def __init__(self, eastings, northings, costs):
        self.eastings = eastings
        self.northings = northings
        self.costs = costs
        self.get_cost = self.create_interpolator()

    def create_interpolator(self):
        points = np.array( (self.eastings.flatten(), self.northings.flatten()) ).T
        values = self.costs.flatten()
        def get_cost(X, Y):
            return griddata( points, values, (X,Y) )
        return get_cost


    def neighbours(self, pt):
        n = []
        # if the last element was a bend (pt.theta > 0) 
        if pt.theta > 0:
            # either another element at that theta increment
            n.append(pt.bend_pt())
            # or a point a minimum straight length away at the same heading
            n.append(pt.straight_pt(min_straight_length))
        else:
            # add all new bend options
            angles = [spacing * 2 * pi / r for r in range(min_bend_radius, max_bend_radius, bend_radius_inc)]
            angles.extend([-1*a for a in angles])
            [n.append(pt.bend_pt(theta)) for theta in angles]
            # and single straight point
            n.append(pt.straight_pt(spacing))

        return filter(self.is_valid, n)

    def is_valid(self, pt):
        return (np.amin(self.eastings) <= pt.x <= np.amax(self.eastings)
                and np.amin(self.northings) <= pt.y <= np.amax(self.northings))

    def cost(self, prev_p, p):
        # this is tricky, we have to use 2D interpolation for each increment of spacing
        # between the 2 pts.
        d = dist(prev_p, p)
        pts = d / spacing + 1
        xs = np.linspace(prev_p.x, p.x, pts)
        ys = np.linspace(prev_p.y, p.y, pts)
        costs = self.get_cost(xs, ys)
        return np.trapz(costs, dx=spacing) / d

def dist(p1, p2):
    return ((p1.x-p2.x)**2 + (p1.y-p2.y)**2)**0.5

def myround(num, div):
   whole = int(num/div)
   if (num%div)/div >= 0.5:
           return div * (whole + 1)
   else:
           return div * whole

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.theta = 0
        self.heading = 0
    
    def bend_pt(self, theta=None):
        if theta is None:
            new_heading = self.heading + self.theta
        else:
            new_heading = self.heading + theta
        new_x = self.x + spacing * sin(new_heading)
        new_y = self.y + spacing * cos(new_heading)
        p = Point(new_x, new_y)
        p.theta = self.theta
        p.heading = new_heading
        return p

    def straight_pt(self, dist):
        new_x = self.x + dist * sin(self.heading)
        new_y = self.y + dist * cos(self.heading)
        p = Point(new_x, new_y)
        p.theta = 0
        p.heading = self.heading
        return p

    def __hash__(self):
        return hash((myround(self.x, spacing), myround(self.y, spacing), myround(self.heading, heading_tol)))

    def __eq__(self, p):
        return self.__hash__() == p.__hash__()

def a_star_pipe(graph, start, goal):
    frontier = PriorityQueue()
    frontier.put(start, 0)
    came_from = {}
    cost_so_far = {}
    came_from[start] = None
    cost_so_far[start] = 0
    
    while not frontier.empty():
        current = frontier.get()
        
        if current == goal:
            break
        if current.x > 30 and current.y > 30:
            break
        
        for next in graph.neighbours(current):
            new_cost = cost_so_far[current] + graph.cost(current, next)
            if next not in cost_so_far or new_cost < cost_so_far[next]:
                cost_so_far[next] = new_cost
                priority = new_cost + heuristic(goal, next)
                frontier.put(next, priority)
                came_from[next] = current
    
    return came_from, cost_so_far

# create grid
X = np.arange(0,50,10)
Y = np.arange(0,50,10)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = (np.sin(R) + 1) * 10
graph = Grid(X, Y, Z)

start = Point(2.4, 6.7)
start.heading = pi/2
goal = Point(44.3, 47.1)
came_from, cost_so_far = a_star_pipe(graph, start, goal)
for p in came_from:
    print(p.x, p.y)