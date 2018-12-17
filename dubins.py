from math import pi, sin, cos, atan, acos, atan2

import numpy as np

from a_star import Point
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

def paths(start, goal, min_bend_radius, min_straight_length, heading_tol, spacing):
    """ Find the shortest path to the goal from a start point
    respecting the minimum bend radius and minimum straight length.
    Points are returned along the path at the spacing requested.

    Paths assessed are increasingly complex until a valid path is
    found. Increased complexity general equates to more bends. Paths
    are assessed in the following order:

        1.  Straight line to goal of minimum straight length;
        2.  Straight-Curve-Straight respecting minimum straight length;
        3.  Straight-Curve-Straight-Curve-Straight
    """
    p = []
    p.extend(S_route(start, goal, min_straight_length, heading_tol, spacing))
    if len(p) > 0:
        # we have a straight path
        return p
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

def S_route(start, goal, min_straight_length, heading_tol, spacing):
    # if they're not even the same heading...
    if not abs(start.heading - goal.heading) < heading_tol:
        # bail out
        return []
    # we check if the start ends up on goal if we move it the distance between them
    d = dist(start, goal)
    # if we're on a bend and not a minimum straight length away
    if start.radius is not None and d < min_straight_length:
        return []
    # no need to add heading since we've already checked that
    p = Point(start.x+d*sin(start.heading), start.y+d*cos(start.heading))
    # if they're not in the same place
    if not p == goal:
        # bail out again
        return []
    return straight_points(start, goal, spacing)

def SCS_route(start, goal, min_straight_length, min_bend_radius, heading_tol, spacing):
    """ THe SCS path is a curve fitted between the start and goal
    points, after both points have been adjusted for the minimum
    straight distance.
    
    For a Straight-Curve-Straight path, the angle between start and
    goal headings must be less than 180 deg or the path will cross over
    itself (not allowed!). This means we need to know the direction to
    turn and the angle .
    """
    
    # we should be using more vectors for this pathfinding stuff
    delta_heading = goal.heading - start.heading
    if delta_heading > pi:
        delta_heading -= 2*pi
    elif delta_heading <= -pi:
        delta_heading += 2*pi
    # adjust the end point for the minimum straight length
    end_tangent_x = goal.x + min_straight_length * sin(goal.heading + pi)
    end_tangent_y = goal.x + min_straight_length * sin(goal.heading + pi)
    end_tangent = Point(end_tangent_x, end_tangent_y, heading=goal.heading)

    return end_tangent.heading

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

def straight_points(start, goal, spacing, leftover=0):
    if leftover == 0:
        s = start
    else:
        leftover_l = spacing - leftover
        s = Point(start.x+leftover_l*sin(start.heading),
                  start.y+leftover_l*cos(start.heading))
    dx = spacing*sin(start.heading)
    dy = spacing*cos(start.heading)
    d = dist(s, goal)
    incs = int(d/spacing)
    xs = np.linspace(start.x, start.x+dx*incs, incs+1)
    ys = np.linspace(start.y, start.y+dy*incs, incs+1)
    return [Point(x, y, heading=start.heading) for x, y in np.stack((xs, ys), axis=-1)]

def dist(p1, p2):
    return ((p1.x-p2.x)**2 + (p1.y-p2.y)**2)**0.5

if __name__ == '__main__':
    
    s = Point(0,0,heading=0)
    g = Point(-10, 10, heading=5*pi/4)

    print(SCS_route(s, g, 2, 1, pi/180, 1))