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
    a = s
    b = Point(s.x+ts*s.dx, s.y+ts*s.dy)
    c = Point(g.x-ts*g.dx, g.y-ts*g.dy)
    d = g
    e = Point(x, y)

    # get distance from vertex to tangent
    tb = v2t(a,b,e,radius)
    te = v2t(b,e,c,radius)
    tc = v2t(e,c,d,radius)

    # get arc lengths
    lab = arc_length(a,b,e,radius)
    lae = arc_length(b,e,c,radius)
    lac = arc_length(e,c,d,radius)

    # we have a number of conditions where we should return infinity
    # if we start at a bend then we must maintain out minimum straight length
    if s.radius != 0 and dist(a,b) - tb < straight:
        return np.inf
    # if not, we must at least get our curve in!
    if s.radius == 0 and dist(a,b) - tb < 0:
        return np.inf
    # if we have a bend at the start, we must also make sure that the next point is the correct direction
    if s.radius != 0 and turn_dir(s, e) * s.radius < 0:
        return np.inf
    # if e lies on line bc the bc must be longer than the minimum 
    if angle(b,e,c) == pi and dist(b,c) - tb - tc < straight:
        return np.inf
    # the straight distance between c and d must be acceptable
    if dist(c, d) - tc < straight:
        return np.inf
    # point e cannot lie on bc outside of b or c
    if angle(b,e,c) == 0:
        return np.inf
    # distance be and ec must be long enough
    if angle(b,e,c) < pi and (dist(b,e)-tb-te < straight or dist(e,c)-te-tc < straight):
        return np.inf
    # if we've got this far we might be ok...
    return dist(a,b) - tb + lab + dist(b,e) - tb - te + lae + dist(e,c) - te - tc + lac + dist(c,d) - tc
    

def route_length(vars, s, g, min_straight, min_bend):
    ts = vars[0]
    tg = vars[1]
    x = vars[2]
    y = vars[3]

    return length(s, g, ts, tg, x, y, min_straight, min_bend)

if __name__ == '__main__':
    s = Point(0, 0, pi/2)
    g = Point(10, 0, -pi/2)
    ts = 10
    tg = 10
    x = 5
    y = 5
    straight = 5
    bend = 2
    vars = np.array([ts, tg, x, y])

    n = minimize(route_length, vars, args=(s, g, straight, bend))
    print(n)




