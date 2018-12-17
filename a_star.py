import heapq
from math import sin, cos, asin, pi, atan2, atan, acos
from copy import deepcopy
import numpy as np
from scipy.interpolate import griddata
from queue import PriorityQueue
import matplotlib.pyplot as plt
import random

spacing = 1
min_straight_length = 10
bend_radius = [2, 4, 6, 8, 10]
theta_tol = 3 * pi / 180

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
        if pt.radius is not None:
            # we keep on this bend or move the minimum straight length
            r = pt.radius
            a = spacing / r
            if pt.direction == 'l':
                h = pt.heading - a
            else:
                h = pt.heading + a
            p = Point(pt.x+spacing*sin(h), pt.y+spacing*cos(h), heading=h, radius=r,
                      direction=pt.direction, previous=pt)
            p.cost = pt.cost + self.cost(pt, p)
            n.append(p)
            p = Point(pt.x+min_straight_length*sin(pt.heading),
                      pt.y+min_straight_length*cos(pt.heading), heading=h, previous=pt)
            p.cost = pt.cost + self.cost(pt, p)
            n.append(p)
        else:
            # we can try every bend angle
            for r in bend_radius:
                a = spacing / r
                for f, d in [(1, 'r'), (-1, 'l')]:
                    h = pt.heading + f * a
                    new_x = pt.x+spacing*sin(h)
                    new_y = pt.y+spacing*cos(h)
                    p = Point(new_x, new_y, previous=pt, heading=h, radius=r)
                    p.cost = pt.cost + self.cost(pt, p)
                    n.append(p)
            new_x = pt.x+spacing*sin(pt.heading)
            new_y = pt.y+spacing*cos(pt.heading)
            p = Point(new_x, new_y, heading=pt.heading, previous=pt)
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
        paths = dubins_paths(p, goal, min(bend_radius), min_straight_length, spacing)
        paths.sort(key=lambda x: x[0])
        for path in paths:
            pts = np.array([(p.x, p.y) for p in path[1]])
            xs = pts[:,0]
            ys = pts[:,1]
            costs = self.get_cost(xs, ys)
            cost = np.trapz(costs, dx=spacing)
            if not np.isnan(cost) and not np.isinf(cost):
                break
        return cost

class Point:

    def __init__(self, x, y, previous=None, cost=None, heading=None, radius=None, direction=None):
        self.x = x
        self.y = y
        self.cost = cost
        self.previous = previous
        self.heading = heading
        self.radius = radius
        self.direction = direction

    def __eq__(self, p):
        if dist(self, p) < spacing:
            if self.heading is not None and p.heading is not None:
                e = abs(self.heading - p.heading) < pi * 1/180
            else:
                e = True
        else:
            e = False
        return e
    
    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"x: {self.x}, y: {self.y}"
class Tangent:

    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2
        if p1.heading == p2.heading:
            self.heading = self.p1.heading
        else:
            raise ValueError('Both points on tangent must be parallel.')

class Circle:

    def __init__(self, centre, radius, direction):
        self.centre = centre
        self.radius = radius
        if direction.lower() in ['l', 'r']:
            self.direction = direction.lower()
        else:
            raise ValueError("direction must be 'l' or 'r'.")
    
    def __repr__(self):
        return f"Circle - centre: {self.centre}, radius: {self.radius}, direction: {self.direction}"

class Arc:

    def __init__(self, circle, start, end):
        self.circle = circle
        self.start = start
        self.end = end
    
    def length(self):
        if self.circle.direction == 'l':
            d_heading = self.start.heading - self.end.heading
            if self.start.heading < self.end.heading:
                d_heading += 2 * pi
        else:
            d_heading = self.end.heading - self.start.heading
            if self.end.heading < self.start.heading:
                d_heading += 2 * pi
        return self.circle.radius * d_heading
    
    def points(self, spacing, leftover=0):
        d_heading = self.end.heading - self.start.heading
        heading_inc = spacing / self.circle.radius
        if self.circle.direction == 'r':
            if d_heading < 0:
                d_heading += 2 * pi
            arc_start_heading = self.start.heading - pi/2
            arc_headings = np.arange(arc_start_heading, arc_start_heading + d_heading, heading_inc)
            pt_headings = arc_headings + pi/2
        else:
            if d_heading > 0:
                d_heading -= 2 * pi
            arc_start_heading = self.start.heading + pi/2
            arc_headings = np.arange(arc_start_heading, arc_start_heading + d_heading, -1*heading_inc)
            pt_headings = arc_headings - pi/2
        xs = self.circle.centre.x + self.circle.radius * np.sin(arc_headings)
        ys = self.circle.centre.y + self.circle.radius * np.cos(arc_headings)
        pts = [Point(x, y, heading=h, radius=self.circle.radius, direction=self.circle.direction) 
               for x, y, h in np.stack((xs, ys, pt_headings), axis=-1)]
        return pts

    def __repr__(self):
        return f"Arc - start: {self.start}, end: {self.end}, centre: {self.circle.centre}"

def dubins_paths(start, goal, min_bend_radius, min_straight_length, spacing):
    """ These will only be true Dubins routes in the cases where the
    start is straight and the tangent length exceeds the minimum straight
    length.
    
    The shortest path to the goal must satisfy our regular route 
    constraints or it will return artificially cheap paths. These are:
        1. Minimum straight length;
        2. Fixed radius bends;
        3. Minimum bends radius.

    If we are currently on a bend then we must either:
        a. Continue at that bend radius to a tangent point;
        b. Continue at the same heading for the minium straight
           length.

    If we are on a straight section then our shortest path must be to
    turn immediately at the minimum bend radius.

    When we are close to the goal i.e. when a tangent length is less
    than the minimum straight length, we must add in a third circle
    the minimum straight length away of the minimum bend radius. When
    the minimum straight length reduces to zero we will return a
    conventional CCC path.
    """
    d = dist(start, goal)
    sl = Point(start.x+d*sin(start.heading), start.y+d*cos(start.heading), heading=start.heading)
    if sl == goal:
        # we can get to the end with a straight line
        return [d, straight_points(start, goal, spacing)]
    routes = CSC_routes(start, goal, min_bend_radius, min_straight_length)
    if len(routes) == 0:
        # We have no valid CSC paths because they all violate our minimum straight length
        # constraint. Now we need to compute CSCSC paths. We have far more options here
        # because we have an extra circle. Potential routes are 'RSRSR', 'RSRSL', 'RSLSR',
        # 'RSLSL', 'LSRSR', 'LSRSL', 'LSLSR' and 'LSLSL'. The direction of the mid circle
        # determines it's location between the start and goal.
        routes = CSCSC_routes(start, goal, min_bend_radius, min_straight_length)
    if len(routes) == 0:
        raise Exception('No valid paths to goal calculated.')
    
    paths = [generate_path(route, spacing) for route in routes]
    return paths

def CSC_routes(start, goal, min_bend_rad, min_straight_length):
    """ Compute all options for Circle-Straight-Circle routes. Each
    route is returned as a list of Arc objects. These Arcs can be
    joined up by their end/start points to generate a path. A list of
    routes is returned from this function.
    """
    routes = []
    if start.radius is not None:
        # we are already on a bend and we know the direction (and we must continue it to
        # maintain our constant bend radius criteria). This limits us to 2 possible paths:
        # 'xSL' and 'xSR'.
        cs = make_circle(start)
        for d in ['l', 'r']:
            cg = make_circle(goal, min_bend_rad, d)
            t = tangent(cs, cg)
            if t is not None and dist(t.p1, t.p2) > min_straight_length:
                routes.append([Arc(cs, start, t.p1), Arc(cg, t.p2, goal)])
    else:
        # we are on a straight section so we should expect up to 4 tangents to be
        # returned. These are the 'RSR', 'LSL', 'RSL' and 'LSR'. When a tangent is 'None'
        # it can be discarded as another CCC path would be shorter anyway.
        for d in [('r', 'r'), ('l', 'l'), ('r', 'l'), ('l', 'r')]:
            cs = make_circle(start, min_bend_rad, d[0])
            cg = make_circle(goal, min_bend_rad, d[1])
            t = tangent(cs, cg)
            if t is not None and dist(t.p1, t.p2) > min_straight_length:
                routes.append([Arc(cs, start, t.p1), Arc(cg, t.p2, goal)])
    return routes
    
def CSCSC_routes(start, goal, min_bend_rad, min_straight_length):
    """ Compute all Circle-Straight-Circle-Straight-Circle routes.
    These should only be required when we have no CSC routes i.e. when
    we are close to the goal and a CSC route would violate the minimum
    straight length constraint.
    """
    routes = []
    CSCSC_paths = [('l', 'r', 'l'), ('l', 'r', 'r'), ('l', 'l', 'l'), ('l', 'l', 'r'),
                   ('r', 'r', 'l'), ('r', 'r', 'r'), ('r', 'l', 'l'), ('r', 'l', 'r')]
    # Find the distances C1 (from the centre of the start circle to the centre of the mid
    # circle) and C2 (from the centre of mid circle to the centre of the goal circle) for
    # the inner and outer tangents the mid and goal circles must always be of minimum bend
    # radius.
    C2_o = min_straight_length
    C2_i = (min_straight_length**2 + (2 * min_bend_rad)**2)**0.5
    # but the start circle radius is based on the radius of the start point
    if start.radius is not None:
        # we're on an arc already so we must maintain the bend radius through the first
        # arc
        C1_o = (min_straight_length**2 - (start.radius - min_bend_rad)**2)**0.5
        C1_i = (min_straight_length**2 + (start.radius + min_bend_rad)**2)**0.5
        # and we can also reduce the path options
        CSCSC_paths = [p for p in CSCSC_paths if p[0] == start.direction]
    else:
        # we are on a straight so the start circle must also be of minimum bend radius
        C1_o = C2_o
        C1_i = C2_i

    for CSCSC in CSCSC_paths:
        # find the distance between circle centres for this case
        if CSCSC[0] == CSCSC[1]:
            C1 = C1_o
        else:
            C1 = C1_i
        if CSCSC[1] == CSCSC[2]:
            C2 = C2_o
        else:
            C2 = C2_i
        # get the centre of the start circle
        if start.radius is not None:
            cs = make_circle(start)
        else:
            cs = make_circle(start, min_bend_rad, CSCSC[0])
        # get the centre of the goal circle:
        cg = make_circle(goal, min_bend_rad, CSCSC[2])
        # get the heading from cs to cg
        if cg.centre.y == cs.centre.y:
            heading_cs_cg = pi / 2
        else:
            heading_cs_cg = atan((cg.centre.x-cs.centre.x)/(cg.centre.y-cs.centre.y)) + pi / 2 # -90 to 270 deg
        d = dist(cs.centre, cg.centre)
        # get the angle between vector cs_cg and vector cs_cm
        if abs(C1**2 + d**2 - C2**2) < abs(2 * C1 * d):
            heading_start_mid = acos((C1**2 + d**2 - C2**2) / (2 * C1 * d))
        else:
            continue
        # get the heading from cs to cm
        if CSCSC[1] == 'r':
            heading_cs_cm = heading_cs_cg - heading_start_mid
        else:
            heading_cs_cm = heading_cs_cg + heading_start_mid
        # make a point at the mid circle centre
        mp = Point(cs.centre.x + C1*sin(heading_cs_cm), cs.centre.y + C1*cos(heading_cs_cm))
        cm = Circle(mp, min_bend_rad, CSCSC[1])
        # get tangents
        t1 = tangent(cs, cm)
        t2 = tangent(cm, cg)
        if t1 is not None and t2 is not None:
            routes.append([Arc(cs, start, t1.p1), Arc(cm, t1.p2, t2.p1), Arc(cg, t2.p2, goal)])
    return routes

def make_circle(pt, radius=None, direction=None):
    if pt.radius is not None:
        r = pt.radius
    else:
        r = radius
    if pt.direction is not None:
        d = pt.direction
    else:
        d = direction
    fd = {'l': -1,
          'r': 1}
    x = pt.x + r * sin(pt.heading + fd[d] * pi/2)
    y = pt.y + r * cos(pt.heading + fd[d] * pi/2)
    return Circle(Point(x, y), r, d)

def tangent(c1, c2):
    """ Returns the Tangent to c1 and c2, travelling in the appropriate directions.
    Tangents are returned as Tangent objects
    Modified version of Java algorithm:
    https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Tangents_between_two_circles

    c1 and c2 are tuples describing the circles and tangent required:
        (cx, cy, rad(>0), d(['l'||'r']))
        cx: float
            x coordinate of circle centre
        cy: float
            y coordinate of circle centre
        rad: float
            circle radius
        d: string ('l' or 'r')
            direction of turning around arc    
    """
    d = dist(c1.centre, c2.centre)
    vx = (c2.centre.x - c1.centre.x) / d
    vy = (c2.centre.y - c1.centre.y) / d

    signs = {('r', 'r'): (1, 1),
             ('r', 'l'): (-1, 1),
             ('l', 'r'): (-1, -1),
             ('l', 'l'): (1, -1)}

    (sign1, sign2) = signs[(c1.direction, c2.direction)]

    c = (c1.radius - sign1 * c2.radius) / d
    if c**2 > 1:
        return None

    h = (max(0.0, 1.0 - c**2))**0.5

    nx = vx * c - sign2 * h * vy
    ny = vy * c + sign2 * h * vx

    x1 = c1.centre.x + c1.radius * nx
    y1 = c1.centre.y + c1.radius * ny
    x2 = c2.centre.x + sign1 * c2.radius * nx
    y2 = c2.centre.y + sign1 * c2.radius * ny
    dx = x2-x1
    dy = y2-y1
    if dy == 0:
        if dx > 0:
            t = pi / 2
        else:
            t = 3 * pi / 2
    else:
        t = pi / 2 - atan2(dy, dx)
    if t < 0:
        t += 2*pi
    return Tangent(Point(x1, y1, heading=t), Point(x2, y2, heading=t))

def generate_path(arcs, spacing):
    leftover = 0
    path_length = 0
    pts = []
    for an in range(len(arcs) -1):
        path_length += arcs[an].length()
        pts.extend(arcs[an].points(spacing, leftover))
        leftover += path_length - (len(pts)-1) * spacing
        path_length += dist(arcs[an].end, arcs[an+1].start)
        pts.extend(straight_points(arcs[an].end, arcs[an+1].start, spacing, leftover))
        leftover += path_length - len(pts) * spacing
    path_length += arcs[an+1].length()
    pts.extend(arcs[an+1].points(spacing, leftover))
    return path_length, pts

def round_down(num, divisor):
    return num - (num%divisor)

def straight_points(start, end, spacing, leftover=0):
    start_x = start.x
    start_y = start.y
    if leftover != 0:
        leftover_l = spacing - leftover
        start_x += leftover_l * sin(start.heading)
        start_y += leftover_l * cos(start.heading)
    dx = round_down(end.x - start_x, spacing)
    dy = round_down(end.y - start_y, spacing)
    p = int(dist(Point(start_x, start_y), end) / spacing)
    xs = np.linspace(start_x, start_x + dx, p)
    ys = np.linspace(start_y, start_y + dy, p)
    pts = [Point(x, y, heading=start.heading) for x, y in np.stack((xs, ys), axis=-1)]
    return pts

def dist(p1, p2):
    return ((p1.x-p2.x)**2 + (p1.y-p2.y)**2)**0.5

def myround(num, div):
   whole = int(num/div)
   if (num%div)/div >= 0.5:
        return div * (whole + 1)
   else:
        return div * whole

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
    X = np.linspace(0, 100, 101)
    Y = np.linspace(0, 100, 101)
    X, Y = np.meshgrid(X, Y)
    Z = np.random.uniform(0,100, size=X.shape)
    graph = Grid(X, Y, Z)

    start_type = random.choice(['line', 'arc'])
    if start_type == 'line':
        start = Point(random.randint(10,90), random.randint(10,90), heading=random.uniform(0, 2*pi))
    else:
        start = Point(random.randint(10,90), random.randint(10,90), heading=random.uniform(0, 2*pi), radius=random.choice(bend_radius), direction=random.choice(['l', 'r']))

    goal = Point(random.randint(10,90), random.randint(10,90), heading=random.uniform(0, 2*pi))

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

