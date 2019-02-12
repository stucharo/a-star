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

def con_dist_ab(vars, s, g, r, min_straight_start, min_straight):
    """ Minimum straight distance from a to b"""
    a, b, c, _ = unpack(s, g, vars)
    return round(dist(a,b) - v2t(a,b,c,r) - min_straight_start, 5)

def con_dist_bc(vars, s, g, r, min_straight):
    """ Minimum straight distance from b to e (if e is not inline)"""
    a, b, c, d = unpack(s, g, vars)
    return dist(b,c) - v2t(a,b,c,r) - v2t(b,c,d,r) - min_straight

def con_dist_cd(vars, s, g, r, min_straight_end, min_straight):
    """ Minimum straight distance from c to d"""
    _, b, c, d = unpack(s, g, vars)
    return dist(c,d) - v2t(b,c,d,r) - min_straight_end

def con_angles(vars, s, g, r, min_straight):
    # constraint to make sure all angles are greater than 0
    a, b, c, d = unpack(s, g, vars)
    return angle(a,b,c) * angle(b,c,d)

def con_turn_dir_b(vars, s, g):
    _, _, c, _ = unpack(s, g, vars)
    return turn_dir(s, c) * s.radius

def get_shortest_route(s, g, min_bend, min_straight, min_straight_end):

    min_length = np.inf
    ips = []

    d, ip = get_length_and_vertices(s, g, min_bend, min_straight, min_straight_end, 'ineq', min_straight)
    if d < min_length:
        min_length = d
        ips = ip
    
    return min_length, ips

def get_length_and_vertices(s, g, min_bend, min_straight, min_straight_end, con, min_straight_start, s_rad=None):

    cons = [{'type': 'ineq', 'fun': con_dist_ab, 'args': (s, g, min_bend, min_straight_start, min_straight)},
            {'type': 'ineq', 'fun': con_dist_bc, 'args': (s, g, min_bend, min_straight)},
            {'type': 'ineq', 'fun': con_dist_cd, 'args': (s, g, min_bend, min_straight_end, min_straight)},
            {'type': 'ineq', 'fun': con_angles, 'args': (s, g, min_bend, min_straight)}]
    
    if s_rad is not None:
        cons.extend({'type': 'ineq', 'fun': con_turn_dir_b, 'args': (s, g)})
    
    start_vars = np.array([min_straight, min_straight_end])

    d = minimize(length, start_vars, args=(s, g, min_bend, min_straight),
                    constraints=cons)

    print(d)

    if d.success: 
        l = d.fun    
        _, b, c, _ = unpack(s, g, d.x)
        ips = [b,c]
    else:
        l = np.inf
        ips = []
        
    return l, ips

def length(vars, s, g, min_bend, min_straight):
    # vars is our minimization variables. these are:
    # vars[0] -> ts
    # vars[1] -> tg
    # vars[2] -> x
    # vars[3] -> y
    # this 
    
    a, b, c, d = unpack(s, g, vars)

    #graph(s, g, [b,c])

    # get distance from vertex to tangent
    tb = v2t(a,b,c,min_bend)
    tc = v2t(b,c,d,min_bend)

    # total distance
    if dist(b, c) - tb - tc > min_straight:
        lab = arc_length(a,b,c,min_bend)
        lac = arc_length(b,c,d,min_bend) 
        return dist(a,b) - tb + lab + dist(b,c) - tb - tc + lac + dist(c,d) - tc
    else:
        return np.inf

def unpack(s, g, vars):
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

if __name__ == '__main__':
    import random
    actual = False
    if actual:
        min_rad = 3
        rads = [-15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        start_rad = 13
        min_straight = 9
        min_straight_end = 4
        sx = 94.98831347446426
        sy = 7.36936759612612
        sh = -0.3962243515382049
        gx = 27.63977577647897
        gy = 74.59474572842517
        gh = -0.7216333115982492
    else:
        max_rad = 15
        min_rad = random.randint(2, max_rad)
        rads = []
        rads.extend(range(-max_rad, -min_rad+1))
        rads.extend([0])
        rads.extend(range(min_rad, max_rad+1))
        start_rad = random.choice(rads)
        sx = random.uniform(0, 100)
        sy = random.uniform(0, 100)
        sh = random.uniform(-pi, pi)
        gx = random.uniform(0, 100)
        gy = random.uniform(0, 100)
        gh = random.uniform(-pi, pi)
        min_straight = random.randint(0, 10)
        min_straight_end = random.randint(0, 10)

    print(f"min_rad = {min_rad}\nrads = {rads}\nstart_rad = {start_rad}\nmin_straight = {min_straight}\nmin_straight_end = {min_straight_end}")
    print(f"sx = {sx}\nsy = {sy}\nsh = {sh}\ngx = {gx}\ngy = {gy}\ngh = {gh}")

    s = Point(sx, sy, sh, radius=start_rad)
    g = Point(gx, gy, gh)

    r = get_shortest_route(s, g, min_rad, min_straight, min_straight_end)

    graph(s, g, r[1])




