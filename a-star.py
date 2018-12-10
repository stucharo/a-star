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
    
    def put(self, item, priority, counter):
        heapq.heappush(self.elements, (priority, counter, item))
    
    def get(self):
        return heapq.heappop(self.elements)[2]

spacing = 1
bend_radius = [1, 2, 3, 4, 5]
heading_tol = pi / 180

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
        # This becomes our heuristic
        # this is tricky, we have to use 2D interpolation for each increment of spacing
        # between the 2 pts.
        d = dist(prev_p, p)
        pts = max(2, myround(d / spacing, spacing) + 1)
        xs = np.linspace(prev_p.x, p.x, pts)
        ys = np.linspace(prev_p.y, p.y, pts)
        costs = self.get_cost(xs, ys)
        cost = np.trapz(costs, dx=d/(pts-1))
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
        return dist(self, p) < spacing and abs(self.heading - p.heading) < heading_tol

def a_star_pipe(graph, start, goal):
    path_ends = PriorityQueue()
    path_ends.put(start, 0, 0)
    counter = 1
    
    while not path_ends.empty():
        # Get path end most likely to be cheapest
        end_point = path_ends.get()
        
        # are we at the end? 
        if end_point == goal:
            break
        
        neighbours = graph.neighbours(end_point)

        for neighbour in neighbours:
            print(neighbour.x, neighbour.y, neighbour.heading)
            priority = neighbour.cost + graph.cost(neighbour, goal)
            path_ends.put(neighbour, priority, counter)
            counter += 1
    
    return end_point

# create grid
X = np.linspace(0, 50, 6)
Y = np.linspace(0, 50, 6)
X, Y = np.meshgrid(X, Y)
Z = np.ones(X.shape)
graph = Grid(X, Y, Z)

start = Point(2.4, 6.7)
start.heading = pi/2
goal = Point(40.4, 6.7)
goal.heading = pi/2

pt = a_star_pipe(graph, start, goal)
while pt.previous is not None:
    print(pt.x, pt.y)
    pt = pt.previous

