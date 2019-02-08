import numpy as np
from scipy.optimize import minimize
from math import pi, sin, cos, asin, acos, atan, atan2, tan

class Point():

    def __init__(self, x, y, heading=0, radius=0):
        self.x = x
        self.y = y
        self.dx = normalize_angle(cos(heading))
        self.dy = normalize_angle(sin(heading))
        self.radius = 0

def normalize_angle(a):
    """ Maintain angle between in range -pi <= a < pi """
    if a > pi:
        return a - 2*pi
    elif a < -pi:
        return a + 2*pi
    else:
        return a

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

def length(s, g, ts, tg, x, y, straight, radius):
    """ Generic length function to calculate path length around max 3 bends
    """
    a = s
    b = Point(s.x+ts*s.dx, s.y+ts*s.dy)
    c = Point(g.x-tg*g.dx, g.y-tg*g.dy)
    d = g
    e = Point(x, y)

    # get distance from vertex to tangent
    tb = v2t(a,b,c,radius)
    tc = v2t(b,c,d,radius)

    # total distance
    if dist(b, c) - tb - tc > straight:
        lab = arc_length(a,b,c,radius)
        lac = arc_length(b,c,d,radius) 
        return dist(a,b) - tb + lab + dist(b,c) - tb - tc + lac + dist(c,d) - tc
    else:
        tb = v2t(a,b,e,radius)
        te = v2t(b,e,c,radius)
        tc = v2t(e,c,d,radius)
        lab = arc_length(a,b,e,radius)
        lae = arc_length(b,e,c,radius)
        lac = arc_length(e,c,d,radius)
        return dist(a,b) - tb + lab + dist(b,e) - tb - te + lae + dist(e,c) - te - tc + lac + dist(c,d) - tc

# The lengths minimization must be constrained to ensure that the rules
# of he routing tool are uphelp

def con_dist_ab(vars, g, s, r, min_straight):
    """ Minimum straight distance from a to b"""
    a = s
    b = Point(s.x+vars[0]*s.dx, s.y+vars[0]*s.dy)
    c = Point(g.x-vars[1]*g.dx, g.y-vars[1]*g.dy)
    d = g
    e = Point(vars[2], vars[3])
    if dist(b,c) - v2t(a,b,c,r) - v2t(b,c,d,r) > min_straight:
        vt = v2t(a,b,c,r)
    else:
        vt = v2t(a,b,e,r)
    return dist(a,b) - vt - min_straight

def con_dist_be(vars, g, s, r, min_straight):
    """ Minimum straight distance from b to e (if e is not inline)"""
    a = s
    b = Point(s.x+vars[0]*s.dx, s.y+vars[0]*s.dy)
    c = Point(g.x-vars[1]*g.dx, g.y-vars[1]*g.dy)
    d = g
    e = Point(vars[2], vars[3])
    if angle(b, e, c) < pi and dist(b,c) - v2t(a,b,c,r) - v2t(b,c,d,r) < min_straight:
        return dist(b, e) - v2t(a,b,e,r) - v2t(b,e,c,r) - min_straight
    else:
        return 1

def con_dist_ec(vars, g, s, r, min_straight):
    """ Minimum straight distance from e to c (if e is not inline)"""
    a = s
    b = Point(s.x+vars[0]*s.dx, s.y+vars[0]*s.dy)
    c = Point(g.x-vars[1]*g.dx, g.y-vars[1]*g.dy)
    d = g
    e = Point(vars[2], vars[3])
    if angle(b, e, c) < pi and dist(b,c) - v2t(a,b,c,r) - v2t(b,c,d,r) < min_straight:
        return dist(e, c) - v2t(b,e,c,r) - v2t(e,c,d,r) - min_straight
    else:
        return 1

def con_dist_cd(vars, g, s, r, min_straight):
    """ Minimum straight distance from c to d"""
    a = s
    b = Point(s.x+vars[0]*s.dx, s.y+vars[0]*s.dy)
    c = Point(g.x-vars[1]*g.dx, g.y-vars[1]*g.dy)
    d = g
    e = Point(vars[2], vars[3])
    if dist(b,c) - v2t(a,b,c,r) - v2t(b,c,d,r) > min_straight:
        vt = v2t(b,c,d,r)
    else:
        vt = v2t(e,c,d,r)
    return dist(c, d) - vt - min_straight

def con_angles(vars, g, s, r, min_straight):
    # constraint to make sure all angles are greater than 0
    a = s
    b = Point(s.x+vars[0]*s.dx, s.y+vars[0]*s.dy)
    c = Point(g.x-vars[1]*g.dx, g.y-vars[1]*g.dy)
    d = g
    e = Point(vars[2], vars[3])
    if dist(b,c) - v2t(a,b,c,r) - v2t(b,c,d,r) > min_straight:
        return angle(a,b,c) * angle(b,c,d)
    else:
        return angle(a,b,e) * angle(b,e,c) * angle(e,c,d)

def route_length(vars, s, g, min_straight, min_bend):
    ts = vars[0]
    tg = vars[1]
    x = vars[2]
    y = vars[3]

    return length(s, g, ts, tg, x, y, min_straight, min_bend)

if __name__ == '__main__':
    s = Point(0, 0, 3*pi/4)
    g = Point(10, 0, pi/4)
    ts = 10
    tg = 10
    x = (s.x+g.x)/2
    y = (s.y+g.y)/2
    straight = 5
    bend = 1
    vars = np.array([ts, tg, x, y])

    cons = [{'type': 'ineq', 'fun': con_angles, 'args': (g, s, bend, straight)},
            {'type': 'ineq', 'fun': con_dist_ab, 'args': (g, s, bend, straight)},
            {'type': 'ineq', 'fun': con_dist_be, 'args': (g, s, bend, straight)},
            {'type': 'ineq', 'fun': con_dist_ec, 'args': (g, s, bend, straight)},
            {'type': 'ineq', 'fun': con_dist_cd, 'args': (g, s, bend, straight)}]

    n = minimize(route_length, vars, args=(s, g, straight, bend), constraints=cons)
    print(n)




