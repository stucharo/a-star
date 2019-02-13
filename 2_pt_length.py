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

def get_shortest_route(s, g, min_bend, min_straight):

    min_length = np.inf
    ips = []

    # d, ip = get_length_and_vertices(s, g, min_bend, min_straight)
    # if d < min_length:
    #     min_length = d
    #     ips = ip

    d, ip = get_length_and_vertices(s, g, min_bend, min_straight, 'eq', 0, 4*min_straight)
    if d < min_length:
        min_length = d
        ips = ip
    
    return min_length, ips

def get_length_and_vertices(s, g, min_bend, min_straight, con, start_straight, min_straight_end, start_radius=None):

    v = [min_straight, min_straight]

    cons = [{'type': con, 'fun': con_len_ab, 'args': (s, g, min_bend, start_straight)},
            {'type': 'ineq', 'fun': con_len_bc, 'args': (s, g, min_bend, min_straight)},
            {'type': 'ineq', 'fun': con_len_cd, 'args': (s, g, min_bend, min_straight_end)}]
    
    d = minimize(length, v, args=(s,g,min_bend, min_straight), constraints=cons)

    if d.success:
        ips = [copy_pt(s,s.dx*d.x[0],s.dy*d.x[0],0), copy_pt(g,g.dx*-d.x[1],g.dy*-d.x[1],0)]
        return d.fun, ips
    else:
        return np.inf, []

def con_len_ab(vars, s, g, r, min_straight):
    a, b, c, _ = unpack(vars, s, g)
    t = v2t(a,b,c,r)
    return vars[0] - t - min_straight

def con_len_bc(vars, s, g, r, min_straight):
    a, b, c, d = unpack(vars, s, g)
    tb = v2t(a,b,c,r)
    tc = v2t(b,c,d,r)
    return dist(b,c) - tb - tc - min_straight

def con_len_cd(vars, s, g, r, min_straight_end):
    _, b, c, d = unpack(vars, s, g)
    t = v2t(b,c,d,r)
    return vars[1] - t - min_straight_end

def length(vars, s, g, min_bend, min_straight):
    # vars is our minimization variables. these are:
    # vars[0] -> ts
    # vars[1] -> tg
    
    a, b, c, d = unpack(vars, s, g)

    # get distance from vertex to tangent
    tb = v2t(a,b,c,min_bend)
    tc = v2t(b,c,d,min_bend)

    lab = arc_length(a,b,c,min_bend)
    lac = arc_length(b,c,d,min_bend) 

    return dist(a,b) - tb + lab + dist(b,c) - tb - tc + lac + dist(c,d) - tc

def unpack(vars, s, g):
    a = s
    b = Point(s.x+vars[0]*s.dx, s.y+vars[0]*s.dy)
    c = Point(g.x-vars[1]*g.dx, g.y-vars[1]*g.dy)
    d = g
    return (a, b, c, d)

def graph(s, g, ips):

    pts = [s, g]
    if len(ips) > 0:
        pts[1:1] = ips
        ps = np.array([(p.x, p.y) for p in pts])
        plt.plot(ps[:,0], ps[:,1])
        plt.axis('equal')
    else:
        cx = (s.x+g.x)/2
        cy = (s.y+g.y)/2
        w = max(abs(g.x-s.x), abs(g.y-s.y))
        plt.axis([cx-w/2-10,cx+w/2+10,cy-w/2-10,cy+w/2+10])
    plt.arrow(s.x, s.y, 5*s.dx, 5*s.dy, head_width=0.5, color='green')
    plt.arrow(g.x, g.y, 5*g.dx, 5*g.dy, head_width=0.5, color='red')
    plt.show()

s = Point(0,0,pi/4,radius=-2)
g = Point(100, 0, pi/4)

graph(s, g, get_shortest_route(s, g, 1, 20)[1])