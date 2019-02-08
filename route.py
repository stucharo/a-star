import numpy as np
from scipy.optimize import minimize
from math import pi, sin, cos, asin, acos, atan, atan2, tan
import matplotlib.pyplot as plt

class Point():

    def __init__(self, x, y, heading=0, radius=0):
        self.x = x
        self.y = y
        self.dx = normalize_angle(cos(heading))
        self.dy = normalize_angle(sin(heading))
        self.radius = radius

    def __str__(self):
        return f"({self.x:.2f}, {self.y:.2f}) + ({self.dx:.2f}, {self.dy:.2f})"

    def __repr__(self):
        return f"({self.x:.2f}, {self.y:.2f}) + ({self.dx:.2f}, {self.dy:.2f})"

def normalize_angle(a):
    """ Maintain angle between in range -pi <= a < pi """
    if a > pi:
        return a - 2*pi
    elif a < -pi:
        return a + 2*pi
    else:
        return a

def dist2intersect(p1, p2):
    dv = np.array([[p1.dx, -p2.dx], [p1.dy, -p2.dy]])
    lv = np.array([p2.x-p1.x, p2.y-p1.y])
    try:
        return np.linalg.solve(dv, lv)
    except:
        return np.array([np.inf,np.inf])
 
def dist(p1, p2):
    return ((p2.x-p1.x)**2+(p2.y-p1.y)**2)**0.5

def angle(p1, p2, p3):
    """ Angle between p1 -> p2 -> p3
    Except clause to catch case when all 3 points are on the same line.
    This returns pi e.g. angle of 180 degress
    """
    try:
        return acos((dist(p2,p3)**2 + dist(p1,p2)**2 - dist(p1,p3)**2)
                    / (2 * dist(p2,p3) * dist(p1,p2)))
    except:
        return pi

def v2t(p1, p2, p3, r):
    """ Distance from vertex `p2` to the tangent of circle radius `r`
    through line p1 -> p2 -> p3.
    """
    a = angle(p1, p2, p3)
    if a == pi:
        return 0
    else:
        return r / tan(a / 2)

def arc_length(p1, p2, p3, r):
    a = angle(p1, p2, p3)
    if a == pi:
        return 0
    else:
        return a * r

def turn_dir(p1, p2):
    cp = np.cross([p1.dx, p1.dy], [p2.x-p1.x, p2.y-p1.y])
    if cp < 0:
        return -1
    else:
        return 1

def copy_pt(p, dx, dy, dt):
    in_angle = atan2(p.dy, p.dx)
    new_heading = normalize_angle(in_angle + dt)
    return Point(p.x+dx, p.y+dy, new_heading, radius=p.radius)

# The lengths minimization must be constrained to ensure that the rules
# of he routing tool are upheld

def con_dist_ab(vars, g, s, start_r, r, min_straight_start, min_straight):
    """ Minimum straight distance from a to b"""
    a, b, c, d, e = unpack(s, g, vars)
    if dist(b,c) - v2t(a,b,c,start_r) - v2t(b,c,d,r) > min_straight:
        vt = v2t(a,b,c,start_r)
    else:
        vt = v2t(a,b,e,start_r)
    return dist(a,b) - vt - min_straight_start

def con_dist_be(vars, g, s, r, min_straight):
    """ Minimum straight distance from b to e (if e is not inline)"""
    a, b, c, d, e = unpack(s, g, vars)
    if angle(b, e, c) < pi and dist(b,c) - v2t(a,b,c,r) - v2t(b,c,d,r) < min_straight:
        return dist(b, e) - v2t(a,b,e,r) - v2t(b,e,c,r) - min_straight
    else:
        return 1

def con_dist_ec(vars, g, s, r, min_straight):
    """ Minimum straight distance from e to c (if e is not inline)"""
    a, b, c, d, e = unpack(s, g, vars)
    if angle(b, e, c) < pi and dist(b,c) - v2t(a,b,c,r) - v2t(b,c,d,r) < min_straight:
        return dist(e, c) - v2t(b,e,c,r) - v2t(e,c,d,r) - min_straight
    else:
        return 1

def con_dist_cd(vars, g, s, r, min_straight_end, min_straight):
    """ Minimum straight distance from c to d"""
    a, b, c, d, e = unpack(s, g, vars)
    if dist(b,c) - v2t(a,b,c,r) - v2t(b,c,d,r) > min_straight:
        vt = v2t(b,c,d,r)
    else:
        vt = v2t(e,c,d,r)
    return dist(c, d) - vt - min_straight_end

def con_angles(vars, g, s, r, min_straight):
    # constraint to make sure all angles are greater than 0
    a, b, c, d, e = unpack(s, g, vars)
    if dist(b,c) - v2t(a,b,c,r) - v2t(b,c,d,r) > min_straight:
        return angle(a,b,c) * angle(b,c,d)
    else:
        return angle(a,b,e) * angle(b,e,c) * angle(e,c,d)

def con_turn_dir_b(vars, g, s, r, min_straight):
    a, b, c, d, e = unpack(s, g, vars)
    if dist(b,c) - v2t(a,b,c,abs(s.radius)) - v2t(b,c,d,r) > min_straight:
        return turn_dir(s, c) * s.radius
    else:
        return turn_dir(s, e) * s.radius

def route_length(vars, s, g, min_straight, min_bend):
    ts = vars[0]
    tg = vars[1]
    x = vars[2]
    y = vars[3]

    return length(s, g, ts, tg, x, y, min_straight, min_bend)

def get_shortest_route(s, g, min_bend, min_straight, min_straight_end):
    """
    Parameters
    ----------
    s: `Point`
        Start point for pathfinder
    g: `Point`
        Goal point for pathfinder
    min_bend: `int` or `float`
        Minimum bend radius for route
    min_straight: `int` or `float`
        Minimum straight length between bends along route
    min_straight_end: `int` or `float`
        Minimum straight length at end of route
    
    Returns
    -------
    Bends
        List of tuples containing bend start and end points along the route

    """
    # we're collecting the lengths of various route options. let's do that in a list of 
    # tuples where each tuple is a list of vertices (excluding s and g) and the total 
    # length e.g. [([p11], 60), ([p21, p22], 45)]
    min_length = np.inf
    ips = []
    # we can minimise a route length, if a valid route is available. we need to know which
    # route function to minimise
    # if the points are too close to fit a mbr in at all then the route returns 0
    # check the distance to intersection between start and end vectors
    d2i = dist2intersect(s, g)
    if 0 < d2i[0] < np.inf:
        # the two lines cross and a triangle is possible
        # lets get a vertex at the intersection point
        ip = Point(s.x+d2i[0]*s.dx, s.y+d2i[0]*s.dy)
        # and the distance between the ip and tp for a minimum radius bend
        d2t = v2t(s, ip, g, min_bend)
        # the minimum required distance depends on whether the start is on a bend
        if s.radius == 0:
            min_d = 0
        else:
            min_d = min_straight
        if d2i[0] < min_d + d2t or -d2i[1] < min_straight_end + d2t:
            # we can't even fit in a minimum bend so return an empty list of bend point
            return min_length, ips
        # if we're here we can make a path with at least 1 bend
        # it's length is...
        min_length = d2i[0]-d2t + arc_length(s,ip,g,min_bend) + -1*d2i[1] + d2t
        # and it's tangent points are
        td = turn_dir(s, g)
        st = copy_pt(s, s.dx*(d2i[0]-d2t), s.dy*(d2i[0]-d2t), td*pi/2)
        gt = copy_pt(g, g.dx*(d2i[1]+d2t), g.dy*(d2i[1]+d2t), td*pi/2)
        ips = [ip]

    # we know there is no option for a single bend path if the two vectors don't cross

    # now we can investigate 2 and 3 bend options

    # can we fit a 2 bend path in?
    # there is 3 options depending on whether the start is on a bend
    # first, lets make initial guesses at where a third point might lie
 
    d, ip = get_length_and_vertices(s, g, min_bend, min_straight, min_straight_end, 0)
    if d < min_length:
        min_length = d
        ips = ip
    
    if s.radius != 0:
        # or move a minimum length length away
        d, ip = get_length_and_vertices(s, g, min_bend, min_straight, min_straight_end)

        if d < min_length:
            min_length = d
            ips = ip
    
    # convert

    return min_length, ips

def get_length_and_vertices(s, g, min_bend, min_straight, min_straight_end, min_straight_start=None):

    cons = [{'type': 'ineq', 'fun': con_angles, 'args': (g, s, min_bend, min_straight)},
            {'type': 'ineq', 'fun': con_dist_be, 'args': (g, s, min_bend, min_straight)},
            {'type': 'ineq', 'fun': con_dist_ec, 'args': (g, s, min_bend, min_straight)},
            {'type': 'ineq', 'fun': con_dist_cd, 'args': (g, s, min_bend, min_straight_end, min_straight)}]

    
    if min_straight_start is None:
        min_straight_start = min_straight

    if s.radius == 0:
        cons.extend([{'type': 'ineq', 'fun': con_dist_ab, 'args': (g, s, min_bend, min_bend, min_straight_start, min_straight)}])
    else:
        cons.extend([{'type': 'eq', 'fun': con_dist_ab, 'args': (g, s, abs(s.radius), min_bend, 0, min_straight)},
                     {'type': 'eq', 'fun': con_turn_dir_b, 'args': (g, s, s.radius, min_straight)}])
    
    start_vars = np.array([2*min_straight_start+1, 2*min_straight_end+1, (s.x+g.x)/2, (s.y+g.y)/2])

    d = minimize(length, start_vars, args=(s, g, min_bend, min_straight),
                    constraints=cons)
                    
    _, b, c, _, e = unpack(s, g, d.x)
    if dist(b, c) > min_straight:
        ips = [b, c]
    else:
        ips=[b, e, c]
    return d.fun, ips 

def length(vars, s, g, min_bend, min_straight):
    # vars is our minimization variables. these are:
    # vars[0] -> ts
    # vars[1] -> tg
    # vars[2] -> x
    # vars[3] -> y
    # this 
    
    a, b, c, d, e = unpack(s, g, vars)

    if s.radius == 0:
        start_radius = min_bend
    else:
        start_radius = abs(s.radius)

    # get distance from vertex to tangent
    tb = v2t(a,b,c,start_radius)
    tc = v2t(b,c,d,min_bend)

    # total distance
    if dist(b, c) - tb - tc > min_straight:
        lab = arc_length(a,b,c,start_radius)
        lac = arc_length(b,c,d,min_bend) 
        return dist(a,b) - tb + lab + dist(b,c) - tb - tc + lac + dist(c,d) - tc
    else:
        tb = v2t(a,b,e,start_radius)
        te = v2t(b,e,c,min_bend)
        tc = v2t(e,c,d,min_bend)
        lab = arc_length(a,b,e,start_radius)
        lae = arc_length(b,e,c,min_bend)
        lac = arc_length(e,c,d,min_bend)
        return dist(a,b) - tb + lab + dist(b,e) - tb - te + lae + dist(e,c) - te - tc + lac + dist(c,d) - tc

def unpack(s, g, vars):
    a = s
    b = Point(s.x+vars[0]*s.dx, s.y+vars[0]*s.dy)
    c = Point(g.x-vars[1]*g.dx, g.y-vars[1]*g.dy)
    d = g
    e = Point(vars[2], vars[3])
    return (a, b, c, d, e)

if __name__ == '__main__':
    s = Point(0, 0, -3*pi/4, radius=10)
    g = Point(10, 0, -pi/4)
    min_bend = 1
    min_straight = 10
    min_straight_end = 10

    r = get_shortest_route(s, g, min_bend, min_straight, min_straight_end)

    pts = [s, g]
    if len(r[1]) > 0:
        pts[1:1] = r[1]
    ps = np.array([(p.x, p.y) for p in pts])
    plt.plot(ps[:,0], ps[:,1])
    plt.axis('equal')
    plt.show()




